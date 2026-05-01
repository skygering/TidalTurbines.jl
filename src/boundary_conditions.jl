export initial_condition_tidal_surge,
       WaveMaker,
       make_initial_condition_tidal_surge,
       make_dirichlet_state,
       make_boundary_condition_flather_left,
       make_boundary_condition_flather_right

# -------------------------------------------------------------------
# Wave maker
# -------------------------------------------------------------------

struct WaveMaker
    A::Float64
    ω::Float64
end

@inline (wm::WaveMaker)(t, eq) = eq.H0 + wm.A * sin(wm.ω * t)

# -------------------------------------------------------------------
# Initial condition factory
# -------------------------------------------------------------------

function make_initial_condition_tidal_surge(bathymetry::AbstractBathymetry)
    return @inline function (x, t, equations::ShallowWaterEquations2D)
        v1 = 0.0
        v2 = 0.0

        x1, x2 = x
        b = bathymetry(x1, x2)

        h = max(equations.threshold_limiter, equations.H0 - b)

        return SVector(h, h * v1, h * v2, b)
    end
end

# -------------------------------------------------------------------
# Boundary condition factories
# -------------------------------------------------------------------

function make_dirichlet_state(η_func, bathymetry)
    return @inline function (x, t, equations::ShallowWaterEquations2D)

        x1, x2 = x
        b = bathymetry(x1, x2)

        η = η_func(t, equations)

        # reconstruct water depth
        h = max(equations.threshold_limiter, η - b)

        # velocity choice
        v1 = 0.0
        v2 = 0.0

        return SVector(h, h*v1, h*v2, b)
    end
end


function make_boundary_condition_flather_left(
    wm::WaveMaker,
    bathymetry::AbstractBathymetry
)
    return @inline function (u_inner, normal_direction, x, t,
                             surface_flux_functions,
                             equations::ShallowWaterEquations2D)

        η_ext = wm(t, equations)

        return boundary_condition_flather(
            η_ext,
            u_inner, normal_direction, x, t,
            surface_flux_functions, equations,
            bathymetry
        )
    end
end


# function make_boundary_condition_flather_right(
#     bathymetry::AbstractBathymetry
# )
#     return @inline function (u_inner, normal_direction, x, t,
#                              surface_flux_functions,
#                              equations::ShallowWaterEquations2D)

#         return boundary_condition_flather(
#             equations.H0,
#             u_inner, normal_direction, x, t,
#             surface_flux_functions, equations,
#             bathymetry
#         )
#     end
# end

# -------------------------------------------------------------------
# Core Flather condition
# -------------------------------------------------------------------

@inline function boundary_condition_flather(η_ext,
    u_inner, normal_direction, x, t,
    surface_flux_functions,
    equations::ShallowWaterEquations2D,
    bathymetry::AbstractBathymetry
)
    sf, ncf = surface_flux_functions

    h, hv1, hv2, b_inner = u_inner
    nx, ny = normal_direction

    # recompute bathymetry at boundary location
    x1, x2 = x
    b = bathymetry(x1, x2)

    hsafe = max(h, 1e-6)

    v1 = hv1 / hsafe
    v2 = hv2 / hsafe
    vn = v1 * nx + v2 * ny

    η_inner = h + b_inner
    h_ext = max(equations.threshold_limiter, η_ext - b)

    c = sqrt(equations.gravity / hsafe)
    vn_ext = vn + c * (η_inner - η_ext)

    tx, ty = -ny, nx
    vt = v1 * tx + v2 * ty

    v1o = vn_ext * nx + vt * tx
    v2o = vn_ext * ny + vt * ty

    u_outer = SVector(h_ext, h_ext * v1o, h_ext * v2o, b)

    return sf(u_inner, u_outer, normal_direction, equations),
           ncf(u_inner, u_outer, normal_direction, equations)
end