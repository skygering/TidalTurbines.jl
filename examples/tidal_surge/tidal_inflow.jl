using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using Trixi2Vtk
using TidalTurbines

rm("examples/tidal_surge/out", recursive=true)
mkpath("examples/tidal_surge/out")


# (1) create mesh
tidal = newProject("flat_tidal", "examples/tidal_surge")
setPolynomialOrder!(tidal, 1)
setMeshFileFormat!(tidal, "ISM-V2")
HOHQMesh.getModelDict(tidal)

D = 40
ymin = 0.0
xmin = 0.0
xmax = 1280.0 / D 
ymax = 480.0 / D
h0 = -50 / D
hshore = -5 / D
bounds = [ymax, ymin, xmin, xmax] # [top, left, bottom, right]
N = [16, 8, 0]
addBackgroundGrid!(tidal, bounds, N)
generate_mesh(tidal)

# (2) equations + ICs + BCs
# bathy = SlopedBathymetry(ymin, ymax, -1.0, -0.1)
headland_rad = 160 / D
bathy = HeadlandBathymetry(xmax / 2, ymax, headland_rad, h0, hshore)


Atide = 0.275
T = 400  # 1 hour in seconds
ω = 2π / T

wm_left  = WaveMaker(Atide, ω)
wm_right = WaveMaker(-Atide, ω)

equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 0.0)

initial_condition = make_initial_condition_tidal_surge(bathy)

boundary_condition = (;
    Bottom = boundary_condition_slip_wall,
    Top    = boundary_condition_slip_wall,
    Left   = make_boundary_condition_flather_left(wm_left, bathy),
    Right  = make_boundary_condition_flather_left(wm_right, bathy),
)

# (3) solver
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
    source_terms = source_terms_manning_friction
)

# (4) time integration
tspan = (0.0, 2 * 3600.0)
ode = semidiscretize(semi, tspan)

# (5) output setput
analysis_callback = AnalysisCallback(semi, interval = 1000)
save_solution = SaveSolutionCallback(
    dt = 0.5,
    output_directory = "examples/tidal_surge/out",
    save_initial_solution = true,
    save_final_solution = true
)
stepsize_callback = StepsizeCallback(cfl = 0.3)
callbacks = CallbackSet(analysis_callback, stepsize_callback, save_solution)

# (6) solve!
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))
sol = solve(ode, SSPRK43(stage_limiter!);
            dt = 0.01,
            ode_default_options()...,
            callback = callbacks,
            adaptive = false)

# (7) transform outputs for paraview           
trixi2vtk("examples/tidal_surge/out/solution_*.h5",
          output_directory = "examples/tidal_surge/out/vtk")