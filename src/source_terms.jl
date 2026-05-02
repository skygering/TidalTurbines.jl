export combine_source_terms, make_sponge_source, source_terms_manning_friction, make_turbine_source

function combine_source_terms(source_terms::Vector)

    return @inline function combined_source(u, x, t, equations)

        S = source_terms[1](u, x, t, equations)

        for i in 2:length(source_terms)
            S += source_terms[i](u, x, t, equations)
        end

        return S
    end
end

function sponge_strength(x, y; Lx, σ_max)

    xnorm = x / Lx

    if xnorm < 0.1
        return σ_max * (1 - xnorm / 0.1)
    elseif xnorm > 0.9
        return σ_max * (1 - (1 - xnorm) / 0.1)
    else
        return 0.0
    end
end

function make_sponge_source(; Lx, σ_max=2.0)

    return @inline function source_sponge(u, x, t, equations::ShallowWaterEquations2D)

        h, hv1, hv2, b = u
        x1, x2 = x

        σ = sponge_strength(x1, x2; Lx=Lx, σ_max=σ_max)

        return SVector(
            zero(eltype(x)),   # mass
            -σ * hv1,          # momentum x
            -σ * hv2,          # momentum y
            zero(eltype(x))    # bathymetry
        )
    end
end

@inline function source_terms_manning_friction(u, x, t, equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    n = 0.001
    h = (h^2 + max(h^2, 1e-8)) / (2 * h)

    Sf = -equations.gravity * n^2 * h^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    return SVector(zero(eltype(x)), Sf * hv_1, Sf * hv_2, zero(eltype(x)))
end

function make_turbine_source(turbines)

    return @inline function source(u, x, t, equations::ShallowWaterEquations2D)

        return source_terms_turbines(u, x, t, equations, turbines)
    end
end

@inline function source_terms_turbines(u, x, t,
                                      equations::ShallowWaterEquations2D,
                                      turbines::Vector{Turbine{T}}) where T
    h, hv1, hv2, _ = u
    zero_T = zero(T)

    h_eff = max(h, equations.threshold_limiter)
    if h_eff <= equations.threshold_limiter
        return SVector(zero_T, zero_T, zero_T, zero_T)
    end

    v1 = hv1 / h_eff
    v2 = hv2 / h_eff
    U  = sqrt(v1^2 + v2^2)

    S_h   = zero_T
    S_hv1 = zero_T
    S_hv2 = zero_T
    S_b   = zero_T

    for turb in turbines

        if h_eff <= turb.h_min
            continue
        end

        dA = turbine_density_single(x, turb)
        if dA <= zero_T
            continue
        end

        At = T(0.25) * T(pi) * turb.D^2
        Ct = Ct_turbine(U, turb)

        Kt = T(0.5) * Ct * At * dA / h_eff

        S_hv1 -= Kt * U * hv1
        S_hv2 -= Kt * U * hv2
    end

    return SVector(S_h, S_hv1, S_hv2, S_b)
end
