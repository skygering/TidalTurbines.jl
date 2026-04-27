export initial_condition_tidal_surge, boundary_condition_tidal_surge

const T_tide = 12.42 * 3600        # M2 tide period (seconds)
const ω = 2π / T_tide              # angular frequency
const A = 1.5                      # meters (adjust per site)

@inline function initial_condition_tidal_surge(x, t, equations::ShallowWaterEquations2D)
    # initially water is at rest
    v1 = 0.0
    v2 = 0.0

    # bottom topography values 
    x1, x2 = x
    b = flat_bathymetry(x1, x2) # TODO: make this modular!

    # shift the water level at dry areas to make sure the water height h stays positive.
    h = max(equations.threshold_limiter, equations.H0 - b)

    # Return the conservative variables
    return SVector(h, h * v1, h * v2, b)
end

@inline function tidal_elevation(t)
    return A * sin(ω * t)
end

@inline function boundary_condition_tidal_surge(u_inner, normal_direction,
                                         x, t, surface_flux_functions,
                                         equations::ShallowWaterEquations2D)

    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # water height from tidal elevation
    b = u_inner[4]
    h_ext = max(equations.threshold_limiter, equations.H0 + tidal_elevation(t) - b)

    # assumtions:
    # (1) shallow water velocity from dispersion relation is v = √(gh)
    # (2) velocity aligned with boundary normal
    v_normal = sqrt(equations.gravity * h_ext)

    u_outer = SVector(h_ext,
                      h_ext * v_normal * normal_direction[1],
                      h_ext * v_normal * normal_direction[2],
                      b)

    flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction, equations)

    return flux, noncons_flux
end