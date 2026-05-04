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
output_directory = "examples/tidal_headland_w_turbines/out"
rm(output_directory, recursive=true, force=true)
mkpath(output_directory)

# problem setup + non_dimensionalization
g = 9.81 # m/s^2
D = 20 # m

# domain values
Lx = 1280 / D
Ly = 480 / D
Lₛ = 0.1 * Lx # sponge layer length
H₀ = -50 / D # bathymetry depth
Hₛ = -5 / D # headland/shelf depth
Rₛ = 160 / D # headland/shelf radius

# tidal values
T = 1 * 60 * 60 * sqrt(g / D) # hr * min/hr * sec/min * sqrt(m/s^2 * 1/m)
ω = 2π / T
Aₜ = 4 * 0.275 / D

# turbine values
u_in = 1 / sqrt(D * g) # cut-in speed
u_rated = 3 / sqrt(D * g) # approx rated velocity
u_out = 5 / sqrt(D * g) # cut-in speed
h_min = 0.02 / D

# sponge layer strength σ [1/s]
# time scale: t ~ L_sponge / u_rated [s]
# non-dimensional control paaram: tσ ~5 is good absorbing layer
σₘ = 5 * u_rated / Lₛ

# solver time values
# T = 0.25 * Lx / Uinf
t_out = round(Int, T / 10)
tspan = (0.0, 2.5 * T)
Δt = 0.5 * sqrt(g / D)

# (1) create mesh
tidal = newProject("tidal_headland_w_turbines", "examples/tidal_headland_w_turbines")
setPolynomialOrder!(tidal, 1)
setMeshFileFormat!(tidal, "ISM-V2")
HOHQMesh.getModelDict(tidal)

bounds = [Ly, 0.0, 0.0, Lx] # [top, left, bottom, right]
N = [32, 16, 0]
addBackgroundGrid!(tidal, bounds, N)
headland_core = newRefinementCenter(
    "headland_core", "smooth",
    [Lx/2, Ly, 0.0],
    Lx/(32*2),
    1.5Rₛ
)
add!(tidal, headland_core)
channel = newRefinementLine("channel", "smooth", [15.0, 7.5, 0.0],
                         [45.0, 7.5, 0.0], Lx/(32*2), 10.0)
add!(tidal, channel)
# fig2 = plotProject!(tidal, GRID + REFINEMENTS)
generate_mesh(tidal)

# (2) equations + ICs + BCs
equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)

# bathy = FlatBathymetry(H₀)
bathy = HeadlandBathymetry(Lx/2, Ly, Rₛ, H₀, Hₛ)

initial_condition = make_initial_condition_tidal_surge(bathy)

wm_left  = WaveMaker(Aₜ, ω)
wm_right = WaveMaker(-Aₜ, ω)
bc_left  = DirichletWaveBC(wm_left, bathy)
bc_right = DirichletWaveBC(wm_right, bathy)
# bc_left  = make_boundary_condition_dirchlet_wave(wm_left, bathy)
# bc_right = make_boundary_condition_dirchlet_wave(wm_right, bathy)

boundary_condition = (;
    Bottom = boundary_condition_slip_wall,
    Top    = boundary_condition_slip_wall,
    Left   = bc_left,
    Right  = bc_right,
)

# (3) make turbines
turbines = [Turbine(; x0 = Lx / 2, y0 = α * Ly, u_rated, u_in, u_out, h_min) for α in LinRange(0.1, 0.4, 3)]

# (4) extra source terms
# friction_source = source_terms_manning_friction
friction_source = FrictionSource()
sponge_source = SpongeSource(Lx, σₘ)
turbine_source = TurbineFriction(turbines)
# sponge_source = make_sponge_source(; Lx, σ_max=σₘ)
# turbine_source = make_turbine_source(turbines)
# source_terms = combine_source_terms((friction_source, sponge_source))
source_terms = CombinedSource((sponge_source, friction_source, turbine_source))
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

mesh_file = joinpath(@__DIR__, "tidal_headland_w_turbines.mesh")
mesh = UnstructuredMesh2D(mesh_file)
semi = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = source_terms
)

# (6) time integration
ode = semidiscretize(semi, tspan)

# (7) output setput
analysis_callback = AnalysisCallback(semi, interval = t_out)

save_solution = SaveSolutionCallback(
    dt = 10 * Δt,
    output_directory = output_directory,
    save_initial_solution = false,
    save_final_solution = true,
)

stepsize_callback = StepsizeCallback(cfl = 0.8)
summary_callback = SummaryCallback()
callbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback, save_solution)

# (8) solve!
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

sol = solve(ode, SSPRK43(stage_limiter!);
            dt = Δt,
            ode_default_options()...,
            callback = callbacks,
            adaptive = false)

summary_callback()

# output for paraview
trixi2vtk(joinpath(output_directory, "solution_*.h5"),
          output_directory = joinpath(output_directory, "vtk"))