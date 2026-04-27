export flat_bathymetry, sloped_bathymetry

flat_bathymetry(x, y) = -25.0

function sloped_bathymetry(x, _)
    x_min = min(x)
    x_max = max(x)
    return -25.0 + 2.0 * (x .- x_min) ./ (x_max .- x_min)
end