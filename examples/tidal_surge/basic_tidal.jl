# this is a flat bottom bathymetry and a sinusoidal M2 tidal period
using TidalTurbines
using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using TrixiBottomTopography
using CairoMakie
using Trixi2Vtk

mkpath("examples/tidal_surge/out")

# (1) generate mesh
tidal_flat = newProject("tidal_flat", "examples/tidal_surge");

setPolynomialOrder!(tidal_flat, 1)
setMeshFileFormat!(tidal_flat, "ISM-V2");
HOHQMesh.getModelDict(tidal_flat);

bounds = [100.0, 0.0, 0.0, 1000.0] 
N = [40, 10, 0]
addBackgroundGrid!(tidal_flat, bounds, N)
generate_mesh(tidal_flat);

# (2) construct solver
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 0.0)

# (2) boundary boundary conditions
initial_condition = initial_condition_tidal_surge;

boundary_condition = (;
    Left   = boundary_condition_wave_maker_m2,
    Right  = boundary_condition_outflow_simple,
    Top    = boundary_condition_slip_wall,
    Bottom = boundary_condition_slip_wall
)

basis = LobattoLegendreBasis(3)

surface_flux = (flux_hll_chen_noelle,
                flux_nonconservative_chen_noelle)

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)


mesh_file = joinpath("examples/tidal_surge", "tidal_flat.mesh")
mesh = UnstructuredMesh2D(mesh_file)

semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    solver;
    boundary_conditions = boundary_condition,
    source_terms = source_terms_manning_friction
)

tspan = (0.0, 0.5 * 12.42 * 3600.0)
ode = semidiscretize(semi, tspan);

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

save_solution = SaveSolutionCallback(
    dt = 100.0,   # 10 minutes
    output_directory = "examples/tidal_surge/out",
    save_initial_solution = true,
    save_final_solution = true
)

stepsize_callback = StepsizeCallback(cfl = 1)

callbacks = CallbackSet(
    analysis_callback,
    stepsize_callback,
    save_solution
)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

sol = solve(ode, SSPRK43(stage_limiter!);
          ode_default_options()...,
          callback = callbacks,
          dt = 1.0,
          adaptive = false)
          

trixi2vtk("examples/tidal_surge/out/solution_*.h5",
          output_directory = "examples/tidal_surge/out/vtk")