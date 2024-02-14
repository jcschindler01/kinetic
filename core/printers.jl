
using Printf

#=
String printers.
=#

## format constants
const printers_COL = 14		## column width
const printers_PRC = 6		## digits of precision (for floats)
const printers_FMT = "f"	## output format key (for floats)

## derived format string for floats
const printers_fmtstring = "%.$(printers_PRC)$(printers_FMT)"


## quickprint
function qp(x::String; col=printers_COL)
	return lpad(x, col-1)*","; 
end

## integers
qp(x::Integer) = qp(string(x))

## floats
eval(:( qp(x::Real) = qp(@sprintf($(printers_fmtstring), x)) ))

## tuple
qp(x::Tuple) = string(qp.(x)...)

## symbol
qp(x::Symbol) = qp(string(x))

# datapoint quickprint labels
function qpl(x::Datapoint)
	s = ""
	for i in 1:(fieldcount(Datapoint)-5)
		s *= qp(fieldnames(Datapoint)[i])
	end
	for k in 1:x.N
		s *= qp("k") * qp("x") * qp("y") * qp("vx") * qp("vy")
	end
	return s
end

# datapoint
function qp(x::Datapoint)
	## outstring
	s = ""
	## parameters and scalars
	for i in 1:(fieldcount(Datapoint)-5)
		s *= qp(getfield(x,i))
	end
	## xy and vxy
	for k in 1:x.N
		s *= qp(k) * qp(x.xy[k,1]) * qp(x.xy[k,2]) * qp(x.vxy[k,1]) * qp(x.vxy[k,2])
	end
	## return
	return s
end


#=
Datapoint to String/File IO.
=#

# datapoint from qp string
function from_qp!(datapoint::Datapoint, s::String)
	## instring
	s = split(s,",")
	## parameters and scalars
	for i in 1:(fieldcount(Datapoint)-5)
		## type and val
		valtype = typeof(getfield(datapoint, i))
		valstring = String(strip(s[i]))
		val = parser(valtype, valstring)
		## assign
		setfield!(datapoint, i, val)
	end
	## xy and vxy
	ss = s[fieldcount(Datapoint)-4:end]
	x  = [parser(Float64, String(strip(ss[5*i-3]))) for i in 1:datapoint.N]
	y  = [parser(Float64, String(strip(ss[5*i-2]))) for i in 1:datapoint.N]
	vx = [parser(Float64, String(strip(ss[5*i-1]))) for i in 1:datapoint.N]
	vy = [parser(Float64, String(strip(ss[5*i-0]))) for i in 1:datapoint.N]
	datapoint.xy  = hcat(x,y)
	datapoint.vxy = hcat(vx,vy)
	## return
	return datapoint
end

from_qp(s::String) = from_qp!(Datapoint(), s)

# datapoint from file
function from_file!(datapoint::Datapoint, filename::String)
	## convert last data string from file
	open(filename, "r") do io
		s = String(chomp(last(readlines(io))))
		from_qp!(datapoint, s)
	end
	## return
	return datapoint
end

from_file(filename::String) = from_file!(Datapoint(), filename)

from_file() = from_file("data.txt")
from_file!(datapoint::Datapoint) = from_file(datapoint, "data.txt")


# datapoint to file
function to_file(datapoint::Datapoint, io::IOStream)
	println(io, qp(datapoint))
end

function to_file(datapoint::Datapoint, filename::String)
	open(filename, "a") do io
		println(io, qp(datapoint))
	end
end

to_file(datapoint::Datapoint) = to_file(datapoint, "data.txt")

## new file
function new_file(datapoint::Datapoint, filename::String)
	open(filename, "w") do io
		println(io, qpl(datapoint))
		println(io,  qp(datapoint))
	end
end

new_file(datapoint::Datapoint) = new_file(datapoint, "data.txt")
