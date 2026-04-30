module TidalTurbines

# Write your package code here.
using Trixi
using TrixiShallowWater
using TrixiBottomTopography

include("bathymetry.jl")
include("boundary_conditions.jl")
include("source_terms.jl")

end
