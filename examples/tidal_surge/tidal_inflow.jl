using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using Trixi2Vtk
using TidalTurbines

mkpath("examples/tidal_surge/out")

# (1) create mesh
tidal = newProject("flat_tidal", "examples/tidal_surge")
setPolynomialOrder!(tidal, 1)
setMeshFileFormat!(tidal, "ISM-V2")
HOHQMesh.getModelDict(tidal)

bounds = [3.402, 0.0, 0.0, 5.488] # [top, left, bottom, right]
N = [16, 8, 0]
addBackgroundGrid!(tidal, bounds, N)
generate_mesh(tidal)

# (2) equations + ICs + BCs
bathy = FlatBathymetry(-1.0)
wm = WaveMaker(0.5, 0.005)
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 0.0)
initial_condition = make_initial_condition_tidal_surge(bathy)
boundary_condition = (;
    Bottom = boundary_condition_slip_wall,
    Top    = boundary_condition_slip_wall,
    Left   = make_boundary_condition_flather_left(wm, bathy),
    Right  = make_boundary_condition_flather_right(bathy),
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
tspan = (0.0, 500)
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