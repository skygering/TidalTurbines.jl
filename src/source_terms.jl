# export combine_source_terms, make_sponge_source, source_terms_manning_friction, make_turbine_source
export CombinedSource, SpongeSource, FrictionSource, TurbineFriction

abstract type AbstractSourceTerm end

struct CombinedSource{S}
    sources::S   # tuple of source structs
end

@inline function (cs::CombinedSource)(u, x, t, equations)
    sources = cs.sources

    S = sources[1](u, x, t, equations)

    @inbounds for i in 2:length(sources)
        S += sources[i](u, x, t, equations)
    end

    return S
end

struct SpongeSource{T} <: AbstractSourceTerm
    Lx::T
    σ_max::T
end

@inline function sponge_strength(x, Lx, σ_max)
    xnorm = x / Lx
    left  = max(0.0, 1 - xnorm / 0.1)
    right = max(0.0, 1 - (1 - xnorm) / 0.1)
    return σ_max * max(left, right)
end

@inline function (s::SpongeSource)(u, x, t, equations)
    h, hv1, hv2, b = u
    x1, x2 = x

    σ = sponge_strength(x1, s.Lx, s.σ_max)

    return SVector(
        zero(eltype(x)),
        -σ * hv1,
        -σ * hv2,
        zero(eltype(x))
    )
end

struct FrictionSource{T} <: AbstractSourceTerm
    n::T
end

FrictionSource() = FrictionSource(0.001)

@inline function (mf::FrictionSource)(u, x, t, equations::ShallowWaterEquations2D)
    h, hv1, hv2, _ = u

    g = equations.gravity
    n2 = mf.n * mf.n
    h2 = h * h

    # wetting correction
    h = (h2 + max(h2, 1e-8)) / (2h)

    vel2 = hv1^2 + hv2^2
    vel  = sqrt(vel2)

    invh = 1 / h
    hpow = invh^(7/3)

    Sf = -g * n2 * hpow * vel

    return SVector(
        zero(eltype(u)),
        Sf * hv1,
        Sf * hv2,
        zero(eltype(u))
    )
end

struct TurbineFriction{T} <: AbstractSourceTerm
    turbines::Vector{Turbine{T}}
end


@inline function (tf::TurbineFriction{T})(u, x, t, equations::ShallowWaterEquations2D) where T
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

    for turb in tf.turbines

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
