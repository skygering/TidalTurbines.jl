# export combine_source_terms, make_sponge_source, source_terms_manning_friction, make_turbine_source
export CombinedSource, SpongeSource, FrictionSource, TurbineFriction

abstract type AbstractSourceTerm end

struct CombinedSource{S<:Tuple}
    sources::S
end

CombinedSource(sources...) = CombinedSource(sources)

@generated function (cs::CombinedSource{S})(u, x, t, equations) where {S}

    N = fieldcount(S)

    terms = [
        :(cs.sources[$i](u, x, t, equations))
        for i in 1:N
    ]

    expr = terms[1]

    for i in 2:N
        expr = :($expr + $(terms[i]))
    end

    return expr
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

@inline function (s::SpongeSource{T})(u, x, t, equations) where T
    h, hv1, hv2, b = u
    x1, x2 = x

    σ = sponge_strength(x1, s.Lx, s.σ_max)
    zeroT = zero(T)
    return SVector(
        zeroT,
        -σ * hv1,
        -σ * hv2,
        zeroT
    )
end

struct FrictionSource{T} <: AbstractSourceTerm
    n::T
end

FrictionSource() = FrictionSource(0.001)

@inline function (mf::FrictionSource{T})(u, x, t, equations::ShallowWaterEquations2D) where T
    h, hv1, hv2, _ = u

    g = equations.gravity
    n2 = mf.n * mf.n

    # wetting correction
    # h2 = h * h
    # h = (h2 + max(h2, 1e-8)) / (2h)

    invh = inv(h)
    hpow = invh * invh * cbrt(invh) # invh^(-7/3)

    vel2 = hv1^2 + hv2^2
    vel  = (vel2)^(1/2)

    Sf = -g * n2 * hpow * vel
    zeroT = zero(T)
    return SVector(
        zeroT,
        Sf * hv1,
        Sf * hv2,
        zeroT
    )
end

struct TurbineFriction{T} <: AbstractSourceTerm
    turbines::Vector{Turbine{T}}
end


@inline function (tf::TurbineFriction{T})(u, x, t,
                                          equations::ShallowWaterEquations2D) where T

    R = promote_type(T, eltype(u))
    zero_R = zero(R)

    h, hv1, hv2, _ = u

    hmin = equations.threshold_limiter
    h_eff = max(h, hmin)

    h_eff <= hmin &&
        return SVector(zero_R, zero_R, zero_R, zero_R)

    invh = inv(h_eff)
    v1 = hv1 * invh
    v2 = hv2 * invh

    U_sq = v1^2 + v2^2
    # # Guard against sqrt(0) — its derivative is infinite, giving NaN under AD
    if U_sq <= zero(U_sq)
        return SVector(zero_R, zero_R, zero_R, zero_R)
    end
    U = (U_sq)^(1/2)

    S_hv1 = zero_R
    S_hv2 = zero_R
    @inbounds for turb in tf.turbines

        h_eff <= turb.h_min && continue

        dA = turbine_density_single(x, turb)
        dA <= zero_R && continue

        Ct = Ct_turbine(U, turb)
        Kt = R(0.5) * Ct * R(turb.At) * dA * invh
        drag = Kt * U
        S_hv1 -= drag * hv1
        S_hv2 -= drag * hv2
    end

    return SVector(zero_R, S_hv1, S_hv2, zero_R)
end