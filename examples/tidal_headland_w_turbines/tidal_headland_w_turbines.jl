using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using Trixi2Vtk
using TidalTurbines

const RESTART_MODE = true

output_directory = "examples/tidal_headland_w_turbines/out"
restart_directory = joinpath(output_directory, "restart")
if !RESTART_MODE
    output_directory = "examples/tidal_headland_w_turbines/out"
    rm(output_directory, recursive=true, force=true)
    mkpath(output_directory)
    restart_directory = joinpath(output_directory, "restart")
    mkpath(restart_directory)   
end

g = 9.81
D = 20

Lx = 1280 / D
Ly = 480 / D
Lₛ = 0.1 * Lx
H₀ = -50 / D
Hₛ = -5 / D
Rₛ = 160 / D

T = 1 * 60 * 60 * sqrt(g / D)
ω = 2π / T
Aₜ = 4 * 0.275 / D

u_in = 1 / sqrt(D * g)
u_rated = 3 / sqrt(D * g)
u_out = 5 / sqrt(D * g)
h_min = 0.02 / D

σₘ = 5 * u_rated / Lₛ

tspan = RESTART_MODE ? (4T, 5T) : (0.0, 4T)
Δt = 0.5 * sqrt(g / D)

tidal = newProject("tidal_headland_w_turbines", "examples/tidal_headland_w_turbines")
setPolynomialOrder!(tidal, 1)
setMeshFileFormat!(tidal, "ISM-V2")
HOHQMesh.getModelDict(tidal)

bounds = [Ly, 0.0, 0.0, Lx]
N = [32, 16, 0]

addBackgroundGrid!(tidal, bounds, N)

headland_core = newRefinementCenter(
    "headland_core", "smooth",
    [Lx/2, Ly, 0.0],
    Lx/(32*2),
    1.5Rₛ
)
add!(tidal, headland_core)

channel = newRefinementLine(
    "channel", "smooth",
    [15.0, 7.5, 0.0],
    [45.0, 7.5, 0.0],
    Lx/(32*2),
    10.0
)
add!(tidal, channel)

generate_mesh(tidal)

equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)

bathy = HeadlandBathymetry(Lx/2, Ly, Rₛ, H₀, Hₛ)

initial_condition = make_initial_condition_tidal_surge(bathy)

wm_left  = WaveMaker(Aₜ, ω)
wm_right = WaveMaker(-Aₜ, ω)

bc_left  = DirichletWaveBC(wm_left, bathy)
bc_right = DirichletWaveBC(wm_right, bathy)

boundary_condition = (;
    Bottom = boundary_condition_slip_wall,
    Top    = boundary_condition_slip_wall,
    Left   = bc_left,
    Right  = bc_right,
)

turbines = [
    Turbine(; x0 = Lx / 2, y0 = α * Ly, u_rated, u_in, u_out, h_min)
    for α in LinRange(0.1, 0.4, 3)
]

friction_source = FrictionSource()
sponge_source = SpongeSource(Lx, σₘ)
turbine_source = TurbineFriction(turbines)

source_terms =  if RESTART_MODE
    CombinedSource((sponge_source, friction_source, turbine_source))
else
     CombinedSource((sponge_source, friction_source, turbine_source))
end

basis = LobattoLegendreBasis(3)

surface_flux = (
    FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
    flux_nonconservative_chen_noelle,
)

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh_file = joinpath(@__DIR__, "tidal_headland_w_turbines.mesh")
mesh = UnstructuredMesh2D(mesh_file)

semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    solver;
    boundary_conditions = boundary_condition,
    source_terms = source_terms
)

save_coeff = RESTART_MODE ? 25 : 25

save_solution = SaveSolutionCallback(
    dt = save_coeff * Δt,
    output_directory = output_directory,
    save_initial_solution = false,
    save_final_solution = true,
)

stepsize_callback = StepsizeCallback(cfl = 0.8)
summary_callback = SummaryCallback()
alive_callback = AliveCallback(alive_interval = 100)

restart_callback = SaveRestartCallback(
    interval = 100000,
    save_final_restart = true,
    output_directory = restart_directory
)

callbacks = if RESTART_MODE
    CallbackSet(summary_callback, stepsize_callback, save_solution, alive_callback)
else
    CallbackSet(summary_callback, stepsize_callback, save_solution, alive_callback, restart_callback)
end

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

if RESTART_MODE
    files = readdir(restart_directory, join=true)
    restart_files = filter(f -> occursin("restart_", basename(f)), files)
    latest_file = maximum(restart_files)

    ode = semidiscretize(semi, tspan)

    integrator = init(
        ode,
        SSPRK43(stage_limiter!);
        dt = Δt,
        callback = callbacks,
        save_everystep = false,
        adaptive = false
    )

    load_timestep!(integrator, latest_file)

    sol = solve!(integrator)

else

    ode = semidiscretize(semi, tspan)

    sol = solve(
        ode,
        SSPRK43(stage_limiter!);
        dt = Δt,
        ode_default_options()...,
        callback = callbacks,
        adaptive = false
    )

end

summary_callback()

if RESTART_MODE
    trixi2vtk(
        joinpath(output_directory, "solution_*.h5"),
        output_directory = joinpath(output_directory, "vtk")
    )
end