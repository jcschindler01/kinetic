
## global params
const COL = 14
const PRC = 6
const FMT = "%._f"

## import
using Printf

## quickprint for strings
function qp(x::String; col=COL)
	return lpad(x, col-1)*","; 
end

## quickprint for integers
qp(x::Integer) = qp(string(x))

## quickprint for floats
#qp(x::Real) = qp(@sprintf("_FMT_", x))       ## (FMT fills _ with PRC)
const qpCMD = replace("""qp(x::Real) = qp(@sprintf("FMT", x))""", "FMT" => (replace(FMT, "_" => string(PRC))))
eval(Meta.parse(qpCMD))
