using HOHQMesh
using OrdinaryDiffEqSSPRK
using CairoMakie
using Trixi2Vtk
using TidalTurbines

function build_base_simulation()
    equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 0.0)

    spline_bathymetry_file = Trixi.download(
        "https://gist.githubusercontent.com/andrewwinters5000/21255c980c4eda5294f91e8dfe6c7e33/raw/1afb73928892774dc3a902e0c46ffd882ef03ee3/monai_bathymetry_data.txt",
        joinpath(@__DIR__, "monai_bathymetry_data.txt")
    )
    bath_spline_struct = BicubicBSpline(spline_bathymetry_file, end_condition = "not-a-knot")
    bathymetry(x, y) = spline_interpolation(bath_spline_struct, x, y)

    @inline function initial_condition_monai_tsunami(x, t, equations::ShallowWaterEquations2D)
        v1 = zero(eltype(x))
        v2 = zero(eltype(x))

        x1, x2 = x
        b = bathymetry(x1, x2)

        h = max(equations.threshold_limiter, equations.H0 - b)
        return SVector(h, h * v1, h * v2, b)
    end

    wavemaker_bc_file = Trixi.download(
        "https://gist.githubusercontent.com/andrewwinters5000/5b11f5f175bddb326d11d8e28398127e/raw/64980e0e4526e0fcd49589b34ee5458b9a1cebff/monai_wavemaker_bc.txt",
        joinpath(@__DIR__, "monai_wavemaker_bc.txt")
    )
    h_spline_struct = CubicBSpline(wavemaker_bc_file; end_condition = "not-a-knot")
    H_from_wave_maker(t) = spline_interpolation(h_spline_struct, t)

    @inline function boundary_condition_wave_maker(u_inner, normal_direction::AbstractVector,
                                                   x, t, surface_flux_functions,
                                                   equations::ShallowWaterEquations2D)
        surface_flux_function, nonconservative_flux_function = surface_flux_functions

        h_ext = max(equations.threshold_limiter, H_from_wave_maker(t) - u_inner[4])

        h0 = 0.13535
        v1_ext = 2 * (sqrt(equations.gravity * h_ext) - sqrt(equations.gravity * h0))

        u_outer = SVector(h_ext, h_ext * v1_ext, zero(eltype(x)), u_inner[4])

        flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction, equations)

        return flux, noncons_flux
    end

    boundary_condition = (
        Bottom = boundary_condition_slip_wall,
        Top    = boundary_condition_slip_wall,
        Right  = boundary_condition_slip_wall,
        Left   = boundary_condition_wave_maker
    )

    volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
    surface_flux = (
        FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                     hydrostatic_reconstruction_chen_noelle),
        flux_nonconservative_chen_noelle
    )

    basis = LobattoLegendreBasis(7)

    indicator_sc = IndicatorHennemannGassnerShallowWater(
        equations, basis,
        alpha_max = 0.5,
        alpha_min = 0.001,
        alpha_smooth = true,
        variable = waterheight
    )

    # volume_integral = VolumeIntegralShockCapturingHG(
    #     indicator_sc;
    #     volume_flux_dg = volume_flux,
    #     volume_flux_fv = surface_flux
    # )
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

    solver = DGSEM(basis, surface_flux, volume_integral)

    mesh_file = joinpath(@__DIR__, "../examples/tsunami/monai_shore.mesh")
    mesh = UnstructuredMesh2D(mesh_file)

    return (; mesh, equations, solver, initial_condition = initial_condition_monai_tsunami,
             boundary_condition)
end

function build_semi(p, base)
    turbines = make_turbines(p)

    local_source_terms = (u, x, t, equations) ->
        source_terms_manning_plus_turbines_param(u, x, t, equations, turbines)

    semi = SemidiscretizationHyperbolic(
        base.mesh,
        base.equations,
        base.initial_condition,
        base.solver;
        boundary_conditions = base.boundary_condition,
        source_terms = local_source_terms
    )

    return semi, turbines
end

function objective(p, base; tspan = (0.0, 22.5), saveat = 0.05, rho = 1000.0)
    semi, turbines = build_semi(p, base)
    ode = semidiscretize(semi, tspan)

    callbacks = CallbackSet()

    power_callback = SimplePowerOutputCallback(semi, turbines;
                                           filename = "out/turbine_power.csv",
                                           dt = 0.05,
                                           rho = 1000.0)

    callbacks = CallbackSet(power_callback);

    stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

    sol = solve(ode, SSPRK43(stage_limiter!); dt = 1.0,
          ode_default_options()..., callback = callbacks, adaptive = false, saveat = saveat, save_everystep = false);
    
    E, _ = compute_total_turbine_energy(sol, semi, turbines; rho = rho)

    return -E
end

base = build_base_simulation()

p0 = [3.80, 1.70]

J0 = objective(p0, base; tspan = (0.0, 1), saveat = 0.05)

println("objective(p0) = ", J0)

using ForwardDiff
g = ForwardDiff.gradient(p -> objective(p, base), p0)
println("gradient at p0 = ", g)