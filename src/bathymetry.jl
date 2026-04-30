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
    x_min::Float64
    x_max::Float64
    b_min::Float64
    b_max::Float64
end

@inline function (b::SlopedBathymetry)(x, y)
    ξ = (x - b.x_min) / (b.x_max - b.x_min)
    return b.b_min + ξ * (b.b_max - b.b_min)
end

# -------------------------
# Spline bathymetry
# -------------------------
struct SplineBathymetry{T} <: AbstractBathymetry
    spline_struct::T
end

@inline (b::SplineBathymetry)(x, y) =
    spline_interpolation(b.spline_struct, x, y)
