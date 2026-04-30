export source_terms_manning_friction

using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using TrixiBottomTopography
using CairoMakie
using Trixi2Vtk


@inline function Cp_from_Ct(Ct)
    Ct_clamped = clamp(Ct, 0.0, 1.0)
    return 0.5 * (1.0 + sqrt(1.0 - Ct_clamped)) * Ct_clamped
end

const Ct_rated_val = 0.516
const Cp_rated_val = Cp_from_Ct(Ct_rated_val)

function make_turbines(p)
    x0, y0 = p

    return [(
        x0 = x0,
        y0 = y0,
        D = 0.30,
        r_support = 0.15,
        Ct_rated = Ct_rated_val,
        Cp_rated = Cp_rated_val,
        u_in = 1.00,
        u_rated = 3.05,
        u_out = 6.00,
        h_min = 0.00
    )]
end

const XI_BUMP = 1.45661

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

@inline function turbine_density_single(x, turb)
    x1, x2 = x
    ψx = bump_1d(x1, turb.x0, turb.r_support)
    ψy = bump_1d(x2, turb.y0, turb.r_support)
    return ψx * ψy / (XI_BUMP * turb.r_support^2)
end

@inline function Ct_turbine(U, turb)
    u_in      = turb.u_in       # 1.0
    u_rated   = turb.u_rated    # 3.05
    u_out     = turb.u_out
    Ct_rated  = turb.Ct_rated   # 0.516

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

@inline function Cp_turbine(U, turb)
    Ct = Ct_turbine(U, turb)
    return Cp_from_Ct(Ct)
end

@inline function Cp_capped_above_rated(U, turb)
    return turb.Cp_rated * (turb.u_rated / U)^3
end

@inline function Ct_from_Cp_cardano(Cp)
    Cp_clamped = clamp(Cp, 0.0, 16.0 / 27.0)
    θ = acos(1.0 - 27.0 * Cp_clamped / 8.0) / 3.0
    a = (2.0 / 3.0) * (1.0 - cos(θ))
    return 4.0 * a * (1.0 - a)
end

# -------------------------------------------------------------------
# Combined source terms:
#   Manning friction + distributed turbine drag
# -------------------------------------------------------------------
# @inline function source_terms_manning_plus_turbines(u, x, t,
#                                                     equations::ShallowWaterEquations2D)
#     h, hv_1, hv_2, _ = u

#     T = eltype(x)
#     zero_T = zero(T)

#     # -------------------------------
#     # Manning friction
#     # -------------------------------
#     n = 0.001
#     h_eff = max(h, equations.threshold_limiter)

#     h_fric = (h_eff^2 + max(h_eff^2, 1e-8)) / (2 * h_eff)
#     Sf = -equations.gravity * n^2 * h_fric^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

#     S_h   = zero_T
#     S_hv1 = Sf * hv_1
#     S_hv2 = Sf * hv_2
#     S_b   = zero_T

#     # -------------------------------
#     # Turbine drag from all turbines
#     # -------------------------------
#     if h_eff > equations.threshold_limiter
#         v1 = hv_1 / h_eff
#         v2 = hv_2 / h_eff
#         U  = sqrt(v1^2 + v2^2)

#         for turb in turbines
#             if h_eff > turb.h_min
#                 dA = turbine_density_single(x, turb)
#                 if dA > 0
#                     At = 0.25 * pi * turb.D^2
#                     Ct = Ct_turbine(U, turb)

#                     Kt = 0.5 * Ct * At * dA / h_eff

#                     S_hv1 -= Kt * U * hv_1
#                     S_hv2 -= Kt * U * hv_2
#                 end
#             end
#         end
#     end

#     return SVector(S_h, S_hv1, S_hv2, S_b)
# end

@inline function source_terms_manning_plus_turbines_param(u, x, t,
                                                          equations::ShallowWaterEquations2D,
                                                          turbines)
    h, hv_1, hv_2, _ = u

    T = eltype(u)
    zero_T = zero(T)

    # Manning friction
    n = 0.001
    h_eff = max(h, equations.threshold_limiter)

    h_fric = (h_eff^2 + max(h_eff^2, T(1e-8))) / (2 * h_eff)
    Sf = -equations.gravity * n^2 * h_fric^(-T(7) / T(3)) * sqrt(hv_1^2 + hv_2^2)

    S_h   = zero_T
    S_hv1 = Sf * hv_1
    S_hv2 = Sf * hv_2
    S_b   = zero_T

    # Turbine drag
    if h_eff > equations.threshold_limiter
        v1 = hv_1 / h_eff
        v2 = hv_2 / h_eff
        U  = sqrt(v1^2 + v2^2)

        for turb in turbines
            if h_eff > turb.h_min
                dA = turbine_density_single(x, turb)
                if dA > zero_T
                    At = T(0.25) * T(pi) * turb.D^2
                    Ct = Ct_turbine(U, turb)
                    Kt = T(0.5) * Ct * At * dA / h_eff

                    S_hv1 -= Kt * U * hv_1
                    S_hv2 -= Kt * U * hv_2
                end
            end
        end
    end

    return SVector(S_h, S_hv1, S_hv2, S_b)
end

# -------------------------------------------------------------------
# Turbine power caculations
# -------------------------------------------------------------------
@inline function power_turbine(U, turb)
    A = π * turb.D^2 / 4
    Cp = Cp_turbine(U, turb)
    return 0.5 * 1025.0 * A * Cp * U^3 / 1e6 # MW
end

# -------------------------------------------------------------------
# Turbine power density diagnostics
# -------------------------------------------------------------------
@inline function turbine_power_density_single(u, x, equations, turb;
                                              rho = 1000.0)
    h, hv_1, hv_2, _ = u

    h_eff = max(h, equations.threshold_limiter)
    if h_eff <= turb.h_min
        return 0.0
    end

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

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_terms_manning_plus_turbines);

