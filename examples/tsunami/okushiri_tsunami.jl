#=
This file confirms that we are able to run documentaiton examples from TrixiShallowWater.jl.
Here we run the example of the Okushiri Tsunami:
(https://trixi-framework.github.io/TrixiShallowWater.jl/stable/tutorials/elixir_shallowwater_monai_tsunami/#Okushiri-Tsunami)
=#

using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using TrixiBottomTopography
using CairoMakie
using Trixi2Vtk

mkpath("examples/tsunami/out")

# (1) plot the raw bathymetry file - this file has columns: x(m), y(m), and z(m)
raw_bathymetry_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/305d203c0409d26075aa2993ff367637/raw/df480a6ff63da1916a19820b060abfea83d40dbf/raw_monai_bathymetry.txt",
                                     joinpath(@__DIR__, "raw_monai_bathymetry.txt"));

file = open(raw_bathymetry_file)
lines = readlines(file)[2:end] # remove header
close(file)

nlines = length(lines)
x = zeros(Float64, nlines)
y = zeros(Float64, nlines)
z = zeros(Float64, nlines)

for (j, line) in enumerate(lines)
    current_line = split(line)
    x[j] = parse(Float64, current_line[1])
    y[j] = parse(Float64, current_line[2])
    z[j] = -parse(Float64, current_line[3])
end

fig1 = Figure()

ax1 = Axis3(fig1[1, 1],
    xlabel = "x [m]",
    ylabel = "y [m]",
    zlabel = "z [m]"
)
CairoMakie.surface!(ax1, x, y, z, colormap = :greenbrownterrain)
display(fig1)
save("figures/tsunami/okushiri_tsunami_bath.png", fig1)

# (2) create a mesh 
monai = newProject("monai_shore", "examples/tsunami");
setPolynomialOrder!(monai, 1)
setMeshFileFormat!(monai, "ISM-V2");
HOHQMesh.getModelDict(monai);

bounds = [3.402, 0.0, 0.0, 5.488] # domain size
N = [8, 4, 0] # coarse background grid
addBackgroundGrid!(monai, bounds, N)

# add refinement regions around the island and the wake/shoreline regions behind the island
island = newRefinementCenter("island", "smooth", [3.36, 1.68, 0.0], 0.1, 0.15)
wake = newRefinementLine("wake", "smooth", [3.75, 1.7, 0.0],
                         [4.75, 1.7, 0.0], 0.15, 0.2)
shoreline_top = newRefinementLine("shoreline", "smooth", [4.816, 3.374, 0.0],
                                  [4.83, 2.366, 0.0], 0.15, 0.168)
shoreline_bottom = newRefinementLine("shoreline", "smooth", [4.97, 2.3, 0.0],
                                     [5.32, 1.4, 0.0], 0.075, 0.22);
add!(monai, island)
add!(monai, wake)
add!(monai, shoreline_top)
add!(monai, shoreline_bottom)
fig2 = plotProject!(monai, GRID + REFINEMENTS)
generate_mesh(monai);

# (3) discritize problem -> construct solver
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 0.0)

# new batheymetry file that is in the correct format for TrixiBottomTopography -> could have reformatted above data too
spline_bathymetry_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/21255c980c4eda5294f91e8dfe6c7e33/raw/1afb73928892774dc3a902e0c46ffd882ef03ee3/monai_bathymetry_data.txt",
                                        joinpath(@__DIR__, "monai_bathymetry_data.txt"));
const bath_spline_struct = BicubicBSpline(spline_bathymetry_file, end_condition = "not-a-knot")
bathymetry(x::Float64, y::Float64) = spline_interpolation(bath_spline_struct, x, y);

@inline function initial_condition_monai_tsunami(x, t, equations::ShallowWaterEquations2D)
    # Initially water is at rest
    v1 = 0.0
    v2 = 0.0

    # Bottom topography values are computed from the bicubic spline created above
    x1, x2 = x
    b = bathymetry(x1, x2)

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard zero due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-13 in double precision, is set in the constructor above
    # for the ShallowWaterEquations and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    h = max(equations.threshold_limiter, equations.H0 - b)

    # Return the conservative variables
    return SVector(h, h * v1, h * v2, b)
end
initial_condition = initial_condition_monai_tsunami;

wavemaker_bc_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/5b11f5f175bddb326d11d8e28398127e/raw/64980e0e4526e0fcd49589b34ee5458b9a1cebff/monai_wavemaker_bc.txt",
                                   joinpath(@__DIR__, "monai_wavemaker_bc.txt"));

const h_spline_struct = CubicBSpline(wavemaker_bc_file; end_condition = "not-a-knot")
H_from_wave_maker(t::Float64) = spline_interpolation(h_spline_struct, t);

@inline function boundary_condition_wave_maker(u_inner, normal_direction::AbstractVector,
                                               x, t, surface_flux_functions,
                                               equations::ShallowWaterEquations2D)
    # Extract the numerical flux functions to compute the conservative and nonconservative
    # pieces of the approximation
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # Compute the water height from the wave maker input file data
    # and then clip to avoid negative water heights and division by zero
    h_ext = max(equations.threshold_limiter, H_from_wave_maker(t) - u_inner[4])

    # Compute the incoming velocity as in Eq. (10) of the paper
    # - S. Vater, N. Beisiegel, and J. Behrens (2019)
    #   A limiter-based well-balanced discontinuous Galerkin method for shallow-water flows
    #   with wetting and drying: Triangular grids
    #   [DOI: 10.1002/fld.4762](https://doi.org/10.1002/fld.4762)
    h0 = 0.13535 # reference incident water height converted to meters
    v1_ext = 2 * (sqrt(equations.gravity * h_ext) - sqrt(equations.gravity * h0))

    # Create the external solution state in the conservative variables
    u_outer = SVector(h_ext, h_ext * v1_ext, zero(eltype(x)), u_inner[4])

    # Calculate the boundary flux and nonconservative contributions
    flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)

    noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction,
                                                 equations)

    return flux, noncons_flux
end;

boundary_condition = (; Bottom = boundary_condition_slip_wall,
                      Top = boundary_condition_slip_wall,
                      Right = boundary_condition_slip_wall,
                      Left = boundary_condition_wave_maker);

@inline function source_terms_manning_friction(u, x, t,
                                               equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    n = 0.001 # friction coefficient
    h = (h^2 + max(h^2, 1e-8)) / (2 * h) # desingularization procedure

    # Compute the common friction term
    Sf = -equations.gravity * n^2 * h^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    return SVector(zero(eltype(x)), Sf * hv_1, Sf * hv_2, zero(eltype(x)))
end;

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(7) # polynomial approximation space with degree 7

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh_file = joinpath(@__DIR__, "monai_shore.mesh")

mesh = UnstructuredMesh2D(mesh_file)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_terms_manning_friction);

tspan = (0.0, 22.5)
ode = semidiscretize(semi, tspan);

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)
save_solution = SaveSolutionCallback(dt = 0.5,
                                    output_directory = "examples/tsunami/out",
                                    save_initial_solution = true,
                                    save_final_solution = true)
stepsize_callback = StepsizeCallback(cfl = 0.6)
callbacks = CallbackSet(analysis_callback,
                        stepsize_callback,
                        save_solution);

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))
sol = solve(ode, SSPRK43(stage_limiter!); dt = 1.0,
          ode_default_options()..., callback = callbacks, adaptive = false);

trixi2vtk("examples/tsunami/out/solution_*.h5", output_directory = "examples/tsunami/out/vtk")