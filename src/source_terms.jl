export source_terms_manning_friction

@inline function source_terms_manning_friction(u, x, t,
                                               equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    n = 0.001 # friction coefficient
    h = (h^2 + max(h^2, 1e-8)) / (2 * h) # desingularization procedure

    # Compute the common friction term
    Sf = -equations.gravity * n^2 * h^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    return SVector(zero(eltype(x)), Sf * hv_1, Sf * hv_2, zero(eltype(x)))
end;