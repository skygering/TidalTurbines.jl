export Turbine, turbine_power_density_single, turbine_power_density_total

struct Turbine{T}
    x0::T
    y0::T
    r_support::T
    Ct_rated::T
    Cp_rated::T
    u_in::T
    u_rated::T
    u_out::T
    h_min::T
    D::T
end

# -------------------------------------------------------------------
# Turbine Constructor
# -------------------------------------------------------------------
function Turbine(; 
    x0, y0,
    u_in, u_rated, u_out,
    r_support = 0.50,
    Ct_rated = 0.516,
    Cp_rated = nothing,
    h_min = 0.0,
    D = 1.0,
)
    if isnothing(Cp_rated)
        Cp_rated = Cp_from_Ct(Ct_rated)
    end
    T = promote_type(typeof(x0), typeof(y0), typeof(r_support),
                     typeof(Ct_rated), typeof(Cp_rated),
                     typeof(u_in), typeof(u_rated),
                     typeof(u_out), typeof(h_min), typeof(D))

    return Turbine{T}(T(x0), T(y0), T(r_support),
                      T(Ct_rated), T(Cp_rated),
                      T(u_in), T(u_rated),
                      T(u_out), T(h_min), T(D))
end

# -------------------------------------------------------------------
# Turbine Control / Power / Thrust Functions
# -------------------------------------------------------------------

@inline function Cp_from_Ct(Ct)
    Ct_clamped = clamp(Ct, 0.0, 1.0)
    return 0.5 * (1.0 + sqrt(1.0 - Ct_clamped)) * Ct_clamped
end

@inline function turbine_density_single(x, turb; XI_BUMP = 1.45661)
    x1, x2 = x
    ψx = bump_1d(x1, turb.x0, turb.r_support)
    ψy = bump_1d(x2, turb.y0, turb.r_support)
    return ψx * ψy / (XI_BUMP * turb.r_support^2)
end

@inline function Ct_turbine(U, turb)
    u_in      = turb.u_in
    u_rated   = turb.u_rated
    u_out     = turb.u_out
    Ct_rated  = turb.Ct_rated

    if U <= 0.0
        return 0.0

    elseif U < u_in
        return Ct_rated * exp(5.0 * (U / u_in - 1.0))

    elseif U <= u_rated
        return Ct_rated

    elseif U < u_out
        Cp = Cp_capped_above_rated(U, turb)
        return Ct_from_Cp_cardano(Cp)
    else
        return 0.0
    end
end

@inline function Ct_from_Cp_cardano(Cp)
    Cp_clamped = clamp(Cp, 0.0, 16.0 / 27.0)
    θ = acos(1.0 - 27.0 * Cp_clamped / 8.0) / 3.0
    a = (2.0 / 3.0) * (1.0 - cos(θ))
    return 4.0 * a * (1.0 - a)
end

@inline function Cp_turbine(U, turb)
    Ct = Ct_turbine(U, turb)
    return Cp_from_Ct(Ct)
end

@inline function Cp_capped_above_rated(U, turb)
    return turb.Cp_rated * (turb.u_rated / U)^3
end

@inline function power_turbine(U, turb)
    A = π * turb.D^2 / 4
    Cp = Cp_turbine(U, turb)
    return 0.5 * 1025.0 * A * Cp * U^3 / 1e6 # MW
end

# -------------------------------------------------------------------
# Turbine power density diagnostics
# -------------------------------------------------------------------
@inline function turbine_power_density_single(u, x, equations, turb;
                                              rho = 1000.0*20^3)
    h, hv_1, hv_2, _ = u

    h_eff = max(h, equations.threshold_limiter)
    # if h_eff <= turb.h_min
    #     return 0.0
    # end

    dA = turbine_density_single(x, turb)
    if dA <= 0
        return 0.0
    end

    v1 = hv_1 / h_eff
    v2 = hv_2 / h_eff
    U  = sqrt(v1^2 + v2^2)

    Cp = Cp_turbine(U, turb)
    At = 0.25 * pi * turb.D^2
    # print("  turbine at (", turb.x0, ", ", turb.y0, "): U = ", U, " m/s, Cp = ", Cp)
    return 0.5 * rho * Cp * At * dA * U^3
end

@inline function turbine_power_density_total(u, x, equations, turbines; rho = 1000.0)
    p_total = zero(eltype(u))
    for turb in turbines
        p_total += turbine_power_density_single(u, x, equations, turb; rho = rho)
    end
    return p_total
end

# -------------------------------------------------------------------
# Smooth step function
# -------------------------------------------------------------------
@inline function smoothstep(s)
    if s <= 0
        return 0.0
    elseif s >= 1
        return 1.0
    else
        return s^2 * (3 - 2s)
    end
end

# -------------------------------------------------------------------
# Bump function for turbine footprint
# -------------------------------------------------------------------
@inline function bump_1d(x, p, r)
    ξ = (x - p) / r
    if abs(ξ) < 1
        return exp(1 - 1 / (1 - ξ^2))
    else
        return zero(x)
    end
end
