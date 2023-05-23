
## module
module Kinetic

## exports
export Datapoint
export evolve!, reversal!
export qp, from_qp!, from_qp
export from_file!, from_file, to_file, new_file
export dist, phasedist, area

## include code
include("datapoints.jl")
include("helpers.jl")
include("printers.jl")
include("evolve.jl")
include("physics.jl")
include("ics.jl")

## end module
end
