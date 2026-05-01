module TidalTurbines

using Trixi
using TrixiShallowWater
using TrixiBottomTopography
using HDF5, WriteVTK

include("bathymetry.jl")
include("boundary_conditions.jl")
include("turbines.jl")
include("source_terms.jl")
include("post_processing.jl")
end
