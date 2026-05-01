using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using Trixi2Vtk
using TidalTurbines

#= 
This setup is based on the paper:

Jordan et al., Combining shallow-water and analytical wake models for tidal array micro-siting, 
Journal of Ocean Engineering and Marine Energy (2022)

https://link.springer.com/article/10.1007/s40722-022-00225-2
=#
output_directory = "examples/tidal_headland/out"
rm(output_directory, recursive=true, force=true)
mkpath(output_directory)

output_directory_w_turbines = "examples/tidal_headland/out_w_turbines"
rm(output_directory_w_turbines, recursive=true, force=true)
mkpath(output_directory_w_turbines)

# problem setup + non_dimensionalization

g = 9.81 # m/s^2
D = 30   # m

# domain values
Lx = 1280 / D
Ly = 480 / D
Lₛ = 0.1 * Lx # sponge layer length
H₀ = -50 / D # bathymetry depth
Hₛ = -5 / D # headland/shelf depth
Rₛ = 160 / D # headland/shelf radius

# tidal values
T = 60 * 60 * sqrt(g / D)
ω = 2π / T
Aₜ = 0.275 / D

# turbine values
u_in = 1 / sqrt(D * g) # cut-in speed
u_rated = 3 / sqrt(D * g) # approx rated velocity
u_out = 5 / sqrt(D * g) # cut-in speed
h_min = 0.02 / D

# sponge layer strength σ [1/s]
# time scale: t ~ L_sponge / u_rated [s]
# non-dimensional control paaram: tσ ~5 is good absorbing layer
σₘ = 1 * u_rated / Lₛ

# solver time values
t_out = round(Int, T / 10)
Δt = 0.5 * sqrt(g / D)
tspan = (0.0, 0.25 * T)

# (1) create mesh
tidal = newProject("tidal_headland", "examples/tidal_headland")
setPolynomialOrder!(tidal, 1)
setMeshFileFormat!(tidal, "ISM-V2")
HOHQMesh.getModelDict(tidal)


bounds = [Ly, 0.0, 0.0, Lx] # [top, left, bottom, right]
N = [16, 8, 0]
addBackgroundGrid!(tidal, bounds, N)
generate_mesh(tidal)

# (2) equations + ICs + BCs
equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)

bathy = HeadlandBathymetry(Lx/2, Ly, Rₛ, H₀, Hₛ)

initial_condition = make_initial_condition_tidal_surge(bathy)

wm_left  = WaveMaker(Aₜ, ω)
wm_right = WaveMaker(-Aₜ, ω)
bc_left  = BoundaryConditionDirichlet(make_dirichlet_state(wm_left, bathy))
bc_right = BoundaryConditionDirichlet(make_dirichlet_state(wm_right, bathy))

boundary_condition = (;
    Bottom = boundary_condition_slip_wall,
    Top    = boundary_condition_slip_wall,
    Left   = bc_left,
    Right  = bc_right,
)

# (3) make turbines
make_turbines(p) = [Turbine(; x0 = p[1], y0 = p[2], u_rated,  u_in, u_out, h_min) for p in p_list]
turbines = make_turbines([(Lx / 2, α * Ly) for α in LinRange(0.175, 0.425, 3)])
# turbines = [Turbine(; x0 = Lx / 2, y0 = α * Ly, u_rated, u_in, u_out, h_min) for α in LinRange(0.175, 0.425, 3)]

# (4) extra source terms
sponge_source = make_sponge_source(; Lx, σ_max=σₘ)
friction_source = source_terms_manning_friction
turbine_source = make_turbine_source(turbines)

source_terms = combine_source_terms([sponge_source, friction_source])
source_terms_w_turbine = combine_source_terms([sponge_source, friction_source, turbine_source])

# (5) solver
surface_flux = (
    FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
    flux_nonconservative_chen_noelle,
)
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassnerShallowWater(
    equations, basis,
    alpha_max = 0.5,
    alpha_min = 0.001,
    alpha_smooth = true,
    variable = waterheight
)

volume_integral = VolumeIntegralShockCapturingHG(
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh_file = joinpath(@__DIR__, "tidal_headland.mesh")
mesh = UnstructuredMesh2D(mesh_file)
semi = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = source_terms
)
semi_w_turbines = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = source_terms_w_turbine
)

# (6) time integration
ode = semidiscretize(semi, tspan)
ode_w_turbines = semidiscretize(semi_w_turbines, tspan)

# (7) output setput
analysis_callback = AnalysisCallback(semi, interval = t_out)
analysis_callback_w_turbines = AnalysisCallback(semi_w_turbines, interval = t_out)

save_solution = SaveSolutionCallback(
    dt = 10 * Δt,
    output_directory = output_directory,
    save_initial_solution = true,
    save_final_solution = true,
)

save_solution_w_turbines = SaveSolutionCallback(
    dt = 10 * Δt,
    output_directory = output_directory_w_turbines,
    save_initial_solution = true,
    save_final_solution = true,
)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(analysis_callback, stepsize_callback, save_solution)
callbacks_w_turbine = CallbackSet(analysis_callback_w_turbines, stepsize_callback, save_solution_w_turbines)

# (8) solve!
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

sol = solve(ode, SSPRK43(stage_limiter!);
            dt = Δt,
            ode_default_options()...,
            callback = callbacks,
            adaptive = false)

sol_w_turbines = solve(ode_w_turbines, SSPRK43(stage_limiter!);
            dt = Δt,
            ode_default_options()...,
            callback = callbacks_w_turbine,
            adaptive = false)
cp(joinpath(output_directory, "mesh.h5"), joinpath(output_directory_w_turbines, "mesh.h5"); force=true)


# (9) re-dimensionalize outputs for paraview 
dim_output_directory = "examples/tidal_headland/dimensional_out"
rm(dim_output_directory, recursive=true, force=true)
mkpath(dim_output_directory)
snapshots = load_all_snapshots(output_directory)
dimensionalized_snapshots = dimensionalize.(snapshots)
export_all_snapshots(dimensionalized_snapshots, dim_output_directory)         

dim_output_directory_w_turbines = "examples/tidal_headland/dimensional_out_w_turbines"
rm(dim_output_directory_w_turbines, recursive=true, force=true)
mkpath(dim_output_directory_w_turbines)
snapshots_w_turbines = load_all_snapshots(output_directory_w_turbines)
dimensionalized_snapshots_w_turbines = dimensionalize.(snapshots_w_turbines)
export_all_snapshots(dimensionalized_snapshots_w_turbines, dim_output_directory_w_turbines) 

dim_output_directory_diff = "examples/tidal_headland/dimensional_out_diffs"
rm(dim_output_directory_diff, recursive=true, force=true)
mkpath(dim_output_directory_diff)
diffs = snapshot_differences(dimensionalized_snapshots, dimensionalized_snapshots_w_turbines)
export_all_snapshots(diffs, dim_output_directory_diff)


trixi2vtk(joinpath(dim_output_directory, "solution_*.h5"),
          output_directory = joinpath(dim_output_directory, "vtk"))

trixi2vtk(joinpath(dim_output_directory_w_turbines, "solution_*.h5"),
          output_directory = joinpath(dim_output_directory_w_turbines, "vtk"))

trixi2vtk(joinpath(dim_output_directory_diff, "solution_*.h5"),
          output_directory = joinpath(dim_output_directory_diff, "vtk"))