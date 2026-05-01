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

rm("examples/tidal_headland/out", recursive=true, force=true)
mkpath("examples/tidal_headland/out")

# problem setup + non_dimensionalization

g = 9.81 # m/s^2
D = 40   # m

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
uᵣ = 3 / sqrt(D * g) # approx rated velocity

# sponge layer strength σ [1/s]
# time scale: t ~ L_sponge / uᵣ [s]
# non-dimensional control paaram: tσ ~5 is good absorbing layer
σₘ = 1 * uᵣ / Lₛ

# solver time values
t_out = round(Int, T / 10)
Δt = 0.5 * sqrt(g / D)
tspan = (0.0, T)

# (1) create mesh
tidal = newProject("flat_tidal", "examples/tidal_headland")
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

# (3) extra source terms
sponge = make_sponge_source(; Lx, σ_max=σₘ)
friction = source_terms_manning_friction
source_terms = combine_source_terms([sponge, friction])

# (4) solver
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

mesh_file = joinpath(@__DIR__, "flat_tidal.mesh")
mesh = UnstructuredMesh2D(mesh_file)
semi = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = source_terms
)

# (5) time integration
ode = semidiscretize(semi, tspan)

# (6) output setput
analysis_callback = AnalysisCallback(semi, interval = t_out)
save_solution = SaveSolutionCallback(
    dt = 10 * Δt,
    output_directory = "examples/tidal_headland/out",
    save_initial_solution = true,
    save_final_solution = true
)
stepsize_callback = StepsizeCallback(cfl = 0.3)
callbacks = CallbackSet(analysis_callback, stepsize_callback, save_solution)

# (6) solve!
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))
sol = solve(ode, SSPRK43(stage_limiter!);
            dt = Δt,
            ode_default_options()...,
            callback = callbacks,
            adaptive = false)

# (7) transform outputs for paraview           
trixi2vtk("examples/tidal_headland/out/solution_*.h5",
          output_directory = "examples/tidal_headland/out/vtk")