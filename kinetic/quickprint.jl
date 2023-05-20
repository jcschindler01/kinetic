
## global params
const COL = 14
const PRC = 6
const FMT = "f"

## import
using Printf

## quickprint for strings
function qp(x::String; col=COL)
	return lpad(x, col-1)*","; 
end

## quickprint for integers
qp(x::Integer) = qp(string(x))

## quickprint for floats
fmtstring = "%.$(PRC)$(FMT)"
## qp(x::Real) = qp(@sprintf("_FMT_", x))
eval(:(qp(x::Real) = qp(@sprintf($(fmtstring), x)) ))

## quickprint of tuple
qp(x::Tuple) = string(qp.(x)...)
