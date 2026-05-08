export Turbine

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
    r_support = 0.15,
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
    return 0.5 * (1.0 + (1.0 - Ct_clamped)^(1/2)) * Ct_clamped
end

@inline function turbine_density_single(x, turb; XI_BUMP = 1.45661)
    x1, x2 = x
    ψx = bump_1d(x1, turb.x0, turb.r_support)
    ψy = bump_1d(x2, turb.y0, turb.r_support)
    return ψx * ψy / (XI_BUMP * turb.r_support^2)
end

# @inline function Ct_turbine(U, turb)
#     u_in      = turb.u_in
#     u_rated   = turb.u_rated
#     u_out     = turb.u_out
#     Ct_rated  = turb.Ct_rated

#     if U <= 0.0
#         return 0.0

#     elseif U < u_in
#         return Ct_rated * exp(5.0 * (U / u_in - 1.0))

#     elseif U <= u_rated
#         return Ct_rated

#     elseif U < u_out
#         Cp = Cp_capped_above_rated(U, turb)
#         return Ct_from_Cp_cardano(Cp)
#     else
#         return 0.0
#     end
# end

# @inline function Ct_from_Cp_cardano(Cp)
#     Cp_clamped = clamp(Cp, 0.0, 16.0 / 27.0)
#     θ = acos(1.0 - 27.0 * Cp_clamped / 8.0) / 3.0
#     a = (2.0 / 3.0) * (1.0 - cos(θ))
#     return 4.0 * a * (1.0 - a)
# end

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
# @inline function turbine_power_density_single(u, x, equations, turb;
#                                               rho = 1000.0)
#     h, hv_1, hv_2, _ = u

#     h_eff = max(h, equations.threshold_limiter)
#     if h_eff <= turb.h_min
#         return 0.0
#     end

#     dA = turbine_density_single(x, turb)
#     if dA <= 0
#         return 0.0
#     end

#     v1 = hv_1 / h_eff
#     v2 = hv_2 / h_eff
#     U  = sqrt(v1^2 + v2^2)

#     Cp = Cp_turbine(U, turb)
#     At = 0.25 * pi * turb.D^2
#     # print("  turbine at (", turb.x0, ", ", turb.y0, "): U = ", U, " m/s, Cp = ", Cp)
#     return 0.5 * rho * Cp * At * dA * U^3
# end

@inline function turbine_power_density_single(u, x, equations, turb;
                                              rho = 1000.0)
    h, hv_1, hv_2, _ = u

    h_eff = max(h, equations.threshold_limiter)
    
    # Use smooth indicator instead of if statement
    # Returns ~1 when h_eff > h_min, ~0 when h_eff < h_min
    h_active = 0.5 * (1.0 + tanh(20.0 * (h_eff - turb.h_min)))

    dA = turbine_density_single(x, turb)
    
    # Use smooth max instead of if dA <= 0
    dA_active = 0.5 * (1.0 + tanh(20.0 * dA))

    v1 = hv_1 / h_eff
    v2 = hv_2 / h_eff
    U  = (v1^2 + v2^2)^(1/2)

    Cp = Cp_turbine(U, turb)
    At = 0.25 * pi * turb.D^2
    
    # Multiply by smooth indicators to smoothly turn off power when conditions aren't met
    power_dens = 0.5 * rho * Cp * At * dA * U^3
    
    return h_active * dA_active * power_dens
end


@inline function turbine_power_density_total(u, x, equations, turbines; rho = 1000.0)
    p_total = zero(eltype(u))
    for turb in turbines
        p_total += turbine_power_density_single(u, x, equations, turb; rho = rho)
    end
    return p_total
end

# # -------------------------------------------------------------------
# # Smooth step function
# # -------------------------------------------------------------------
# @inline function smoothstep(s)
#     if s <= 0
#         return 0.0
#     elseif s >= 1
#         return 1.0
#     else
#         return s^2 * (3 - 2s)
#     end
# end

# # -------------------------------------------------------------------
# # Bump function for turbine footprint
# # -------------------------------------------------------------------
# @inline function bump_1d(x, p, r)
#     ξ = (x - p) / r
#     if abs(ξ) < 1
#         return exp(1 - 1 / (1 - ξ^2))
#     else
#         return zero(x)
#     end
# end

@inline function Ct_turbine(U, turb)
    u_in      = turb.u_in
    u_rated   = turb.u_rated
    u_out     = turb.u_out
    Ct_rated  = turb.Ct_rated

    # Region 1: Below cut-in (U < u_in)
    U_ratio_in = U / u_in
    Ct_ramp = Ct_rated * exp(5.0 * (U_ratio_in - 1.0))
    
    # Region 2: Rated region (u_in <= U <= u_rated)
    Ct_rated_region = Ct_rated
    
    # Region 3: Power-controlled region (u_rated < U < u_out)
    Cp = Cp_capped_above_rated(U, turb)
    Cp_safe = clamp(Cp, 0.0, 16.0/27.0)
    arg = clamp(1.0 - 27.0 * Cp_safe / 8.0, -1.0, 1.0)
    θ = acos(arg) / 3.0
    a_temp = (2.0 / 3.0) * (1.0 - cos(θ))
    a_temp = clamp(a_temp, 0.0, 1.0)
    Ct_power_control = 4.0 * a_temp * (1.0 - a_temp)
    
    # Smooth transitions using tanh
    w_in = 0.5 * (1.0 + tanh(20.0 * (U - u_in)))
    Ct_low = Ct_ramp * (1.0 - w_in) + Ct_rated_region * w_in
    
    w_rated = 0.5 * (1.0 + tanh(20.0 * (U - u_rated)))
    Ct_mid = Ct_low * (1.0 - w_rated) + Ct_power_control * w_rated
    
    w_out = 0.5 * (1.0 + tanh(20.0 * (U - u_out)))
    Ct = Ct_mid * (1.0 - w_out)
    
    return Ct
end

@inline function Ct_from_Cp_cardano(Cp)
    # Clamp Cp first, but more carefully
    Cp_min = 0.0
    Cp_max = 16.0 / 27.0
    
    # Use smooth min/max to avoid discontinuous derivatives
    Cp_safe = Cp_max * 2.0 / π * atan(π / 2.0 * Cp / Cp_max)
    
    # Now safe to compute
    arg = 1.0 - 27.0 * Cp_safe / 8.0
    
    # Avoid exact boundaries of acos domain
    arg = clamp(arg, -0.9999, 0.9999)
    
    θ = acos(arg) / 3.0
    a = (2.0 / 3.0) * (1.0 - cos(θ))
    
    return 4.0 * a * (1.0 - a)
end

@inline function bump_1d(x, p, r)
    ξ = (x - p) / r
    ξ_sq = ξ^2
    
    # Use smooth indicator with tanh instead of if statement
    # This smoothly transitions from 1 inside to 0 outside
    indicator = 0.5 * (1.0 + tanh(10.0 * (1.0 - ξ_sq)))
    
    # Main bump function with smooth cutoff
    bump = exp(1.0 - 1.0 / (max(1.0 - ξ_sq, 1e-15))) * indicator
    
    return bump
end