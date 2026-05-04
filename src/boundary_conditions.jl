export initial_condition_tidal_surge,
    make_initial_condition_tidal_surge,
    make_boundary_condition_dirchlet_velocity,
    make_boundary_condition_dirchlet_elevation,
    make_boundary_condition_dirchlet_wave,
    WaveMaker,
    DirichletWaveBC,
    make_dirichlet_state,
    make_boundary_condition_flather_left,
    make_boundary_condition_flather_right

# -------------------------------------------------------------------
# Initial condition
# -------------------------------------------------------------------

function make_initial_condition_tidal_surge(bathymetry::AbstractBathymetry; v1 = 0.0, v2 = 0.0)
    return @inline function (x, t, equations::ShallowWaterEquations2D)
        x1, x2 = x
        b = bathymetry(x1, x2)
        h = max(equations.threshold_limiter, equations.H0 - b)
        return SVector(h, h * v1, h * v2, b)
    end
end

# -------------------------------------------------------------------
# Uniform inflow / flat outflow BCs from:
# Jordan et al., Combining shallow-water and analytical wake models for tidal array micro-siting, 
# Journal of Ocean Engineering and Marine Energy (2022)
# -------------------------------------------------------------------

function make_boundary_condition_dirchlet_velocity(
    u_in::Float64,
    bathymetry::AbstractBathymetry
)
    return @inline function (u_inner, normal_direction, x, t,
                             surface_flux_functions,
                             equations::ShallowWaterEquations2D,
    )
        h, hv1, hv2, b_inner = u_inner
        # get bathymetry
        x1, x2 = x
        b = bathymetry(x1, x2)
        # get elevation
        h = max(equations.threshold_limiter, h)
        # set velocity
        v1 = u_in
        v2 = 0.0
        # get fluxes
        u_outer = SVector(h, h*v1, h*v2, b)
        sf, ncf = surface_flux_functions
        return sf(u_inner, u_outer, normal_direction, equations),
               ncf(u_inner, u_outer, normal_direction, equations)
    end
end

function make_boundary_condition_dirchlet_elevation(
    η_out::Float64,
    bathymetry::AbstractBathymetry
)
    return @inline function (u_inner, normal_direction, x, t,
                             surface_flux_functions,
                             equations::ShallowWaterEquations2D)
        h, hv1, hv2, b_inner = u_inner
        # get bathymetry
        x1, x2 = x
        b = bathymetry(x1, x2)
        # enforce elevation
        h_ext = max(equations.threshold_limiter, η_out - b)
        # KEEP interior velocity
        hsafe = max(h, 1e-6)
        v1 = hv1 / hsafe
        v2 = hv2 / hsafe
        # output
        u_outer = SVector(h_ext, h_ext * v1, h_ext * v2, b)
        sf, ncf = surface_flux_functions
        return sf(u_inner, u_outer, normal_direction, equations),
               ncf(u_inner, u_outer, normal_direction, equations)
    end
end

# -------------------------------------------------------------------
# Tidal inflow / opposing height boundary conditions:
# Jordan et al., Combining shallow-water and analytical wake models for tidal array micro-siting, 
# Journal of Ocean Engineering and Marine Energy (2022)
# -------------------------------------------------------------------

struct WaveMaker
    A::Float64
    ω::Float64
end

@inline (wm::WaveMaker)(t, eq) = eq.H0 + wm.A * sin(wm.ω * t)

# -------------------------------------------------------------------
# Boundary condition factories
# -------------------------------------------------------------------

function make_dirichlet_state(η_func, bathymetry)
    return @inline function (x, t, equations::ShallowWaterEquations2D)

        x1, x2 = x
        b = bathymetry(x1, x2)
        η = η_func(t, equations)
        # reconstruct water depth
        # println("Time: ", t)
        # println(equations.threshold_limiter)
        # println(η)
        # println(b)
        # sleep(0.1)
        # h = η - b
        h = max(equations.threshold_limiter, η - b)
        # println(h)
        # sleep(0.1)


        # velocity choice
        v1 = 0.0
        v2 = 0.0

        return SVector(h, h*v1, h*v2, b)
    end
end

struct DirichletWaveBC{WM, B}
    η_func::WM
    bathy::B
end

@inline function (bc::DirichletWaveBC)(
    u_inner, normal_direction, x, t,
    surface_flux_functions::Tuple{SF, NCF},
    equations::ShallowWaterEquations2D
) where {SF, NCF}
    η_func = bc.η_func
    bathymetry = bc.bathy
    h, hv1, hv2, b_inner = u_inner
    # get bathymetry
    x1, x2 = x
    b = bathymetry(x1, x2)
    # enforce elevation
    η = η_func(t, equations)
    h_ext = max(equations.threshold_limiter, η - b)
    # KEEP interior velocity
    hsafe = max(h, 1e-6)
    v1 = hv1 / hsafe
    v2 = hv2 / hsafe
    # output
    u_outer = SVector(h_ext, h_ext * v1, h_ext * v2, b)
    sf, ncf = surface_flux_functions
    return sf(u_inner, u_outer, normal_direction, equations),
            ncf(u_inner, u_outer, normal_direction, equations)
end


function make_boundary_condition_dirchlet_wave(
    η_func::WaveMaker,
    bathymetry::AbstractBathymetry
)
    return @inline function (u_inner, normal_direction, x, t,
                            surface_flux_functions::Tuple{SF, NCF},
                            equations::ShallowWaterEquations2D
    ) where {SF, NCF}
        h, hv1, hv2, b_inner = u_inner
        # get bathymetry
        x1, x2 = x
        b = bathymetry(x1, x2)
        # enforce elevation
        η = η_func(t, equations)
        h_ext = max(equations.threshold_limiter, η - b)
        # KEEP interior velocity
        hsafe = max(h, 1e-6)
        v1 = hv1 / hsafe
        v2 = hv2 / hsafe
        # output
        u_outer = SVector(h_ext, h_ext * v1, h_ext * v2, b)
        sf, ncf = surface_flux_functions
        return sf(u_inner, u_outer, normal_direction, equations),
               ncf(u_inner, u_outer, normal_direction, equations)
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


function make_boundary_condition_flather_right(
    bathymetry::AbstractBathymetry
)
    return @inline function (u_inner, normal_direction, x, t,
                             surface_flux_functions,
                             equations::ShallowWaterEquations2D)

        return boundary_condition_flather(
            equations.H0,
            u_inner, normal_direction, x, t,
            surface_flux_functions, equations,
            bathymetry
        )
    end
end

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

    # get velocities
    hsafe = max(h, 1e-6)
    v1 = hv1 / hsafe
    v2 = hv2 / hsafe

    # calculate bc
    η_inner = h + b_inner
    h_ext = max(equations.threshold_limiter, η_ext - b)
    vn_ext = 0 + sqrt(equations.gravity / hsafe) * (η_inner - η_ext)

    # reconstruct velocity
    tx, ty = -ny, nx
    vt = v1 * tx + v2 * ty
    v1b = vn_ext * nx + vt * tx
    v2b = vn_ext * ny + vt * ty

    u_outer = SVector(h_ext, h_ext * v1b, h_ext * v2b, b)

    return sf(u_inner, u_outer, normal_direction, equations),
           ncf(u_inner, u_outer, normal_direction, equations)
end