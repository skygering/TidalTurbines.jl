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

output_directory = "examples/tidal_headland_objective/out"
rm(output_directory, recursive=true, force=true)
mkpath(output_directory)

# problem setup + non_dimensionalization
g = 9.81 # m/s^2
D = 20   # m

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
tidal = newProject("tidal_headland", "examples/tidal_headland_objective")
setPolynomialOrder!(tidal, 1)
setMeshFileFormat!(tidal, "ISM-V2")
HOHQMesh.getModelDict(tidal)

bounds = [Ly, 0.0, 0.0, Lx] # [top, left, bottom, right]
N = [16, 8, 0]
addBackgroundGrid!(tidal, bounds, N)
generate_mesh(tidal)

# function that returns solver setup
function build_base_simulation()
    # mesh
    mesh_file = joinpath("examples/tidal_headland_objective", "tidal_headland.mesh")
    mesh = UnstructuredMesh2D(mesh_file)
    
    # solver
    basis = LobattoLegendreBasis(3)
    surface_flux = (
    FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
        flux_nonconservative_chen_noelle,
    )
    volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)
    
    # equations
    equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)
    
    # intitial conditions
    bathy = HeadlandBathymetry(Lx/2, Ly, Rₛ, H₀, Hₛ)
    initial_condition = make_initial_condition_tidal_surge(bathy)
    
    # boundary conditions
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
    return (; mesh, solver, equations, initial_condition, boundary_condition)
end

function build_semi(p_list, base)
    # turbines
    turbines = [Turbine(; x0 = p[1], y0 = p[2], u_rated, u_in, u_out, h_min) for p in p_list]

    # source terms
    sponge_source = make_sponge_source(; Lx, σ_max=σₘ)
    friction_source = source_terms_manning_friction
    turbine_source = make_turbine_source(turbines)
    all_source_terms = combine_source_terms([sponge_source, friction_source, turbine_source])

    semi = SemidiscretizationHyperbolic(
        base.mesh, base.equations, base.initial_condition, base.solver;
        boundary_conditions = base.boundary_condition,
        source_terms = all_source_terms
    )
    return semi, turbines
end

function objective(p, base; tspan, saveat = 0.05, rho = 1000.0)
    semi, turbines = build_semi(p, base)
    ode = semidiscretize(semi, tspan)

    callbacks = CallbackSet()

    stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

    sol = solve(ode, SSPRK43(stage_limiter!);
                dt = Δt,
                ode_default_options()...,
                callback = callbacks,
                adaptive = false,
                saveat)

    E, _ = compute_total_turbine_energy(sol, semi, turbines; rho)
    return -E
end

base = build_base_simulation()
p = [(Lx / 2, 0.25 * Ly)]
T = 2

J0 = objective(p, base; tspan = (0.0, T))

## Optimization challenge
using ForwardDiff
gradient = ForwardDiff.gradient(objective, )