
## module
module Kinetic

## exports
export Datapoint
export evolve!, reversal!
export qp, from_qp!, from_qp
export from_file!, from_file, to_file, new_file
export dist, phasedist, area
export Stau, spatial_macrostate, spatial_log2_qf, S_spatial, S_velocity 
export S_localE, S_vvector

## include code
include("datapoints.jl")
include("helpers.jl")
include("printers.jl")
include("evolve.jl")
include("physics.jl")
include("ics.jl")
include("entropic.jl")

## end module
end

