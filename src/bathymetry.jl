export AbstractBathymetry,
       FlatBathymetry,
       SlopedBathymetry,
       SplineBathymetry

# -------------------------------------------------------------------
# Bathymetry hierarchy
# -------------------------------------------------------------------

abstract type AbstractBathymetry end

# -------------------------
# Flat bathymetry
# -------------------------
struct FlatBathymetry <: AbstractBathymetry
    b0::Float64
end

@inline (b::FlatBathymetry)(x, y) = b.b0

# -------------------------
# Sloped bathymetry
# -------------------------
struct SlopedBathymetry <: AbstractBathymetry
    y_min::Float64
    y_max::Float64
    b_min::Float64
    b_max::Float64
end

@inline function (b::SlopedBathymetry)(x, y)
    ξ = (y - b.y_min) / (b.y_max - b.y_min)
    shape = 0.5 * (1 - cos(2π * ξ))  # 0 at edges, 1 at center
    return b.b_max - shape * (b.b_max - b.b_min)
end

# -------------------------
# Spline bathymetry
# -------------------------
struct SplineBathymetry{T} <: AbstractBathymetry
    spline_struct::T
end

@inline (b::SplineBathymetry)(x, y) =
    spline_interpolation(b.spline_struct, x, y)
