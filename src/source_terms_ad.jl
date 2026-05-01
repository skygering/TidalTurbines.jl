export source_terms_manning_friction

using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using TrixiBottomTopography
using CairoMakie
using Trixi2Vtk

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

# -------------------------------------------------------------------
# Combined source terms:
#   Manning friction + distributed turbine drag
# -------------------------------------------------------------------
@inline function source_terms_manning_plus_turbines(u, x, t,
                                                    equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    T = eltype(x)
    zero_T = zero(T)

    # -------------------------------
    # Manning friction
    # -------------------------------
    n = 0.001
    h_eff = max(h, equations.threshold_limiter)

    h_fric = (h_eff^2 + max(h_eff^2, 1e-8)) / (2 * h_eff)
    Sf = -equations.gravity * n^2 * h_fric^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    S_h   = zero_T
    S_hv1 = Sf * hv_1
    S_hv2 = Sf * hv_2
    S_b   = zero_T

    # -------------------------------
    # Turbine drag from all turbines
    # -------------------------------
    if h_eff > equations.threshold_limiter
        v1 = hv_1 / h_eff
        v2 = hv_2 / h_eff
        U  = sqrt(v1^2 + v2^2)

        for turb in turbines
            if h_eff > turb.h_min
                dA = turbine_density_single(x, turb)
                if dA > 0
                    At = 0.25 * pi * turb.D^2
                    Ct = Ct_turbine(U, turb)

                    Kt = 0.5 * Ct * At * dA / h_eff

                    S_hv1 -= Kt * U * hv_1
                    S_hv2 -= Kt * U * hv_2
                end
            end
        end
    end

    return SVector(S_h, S_hv1, S_hv2, S_b)
end

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
