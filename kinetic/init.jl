
# Usage: julia init.jl -N 50 -i corner

include("./quickprint.jl")
include("./ics.jl")

function help()
	helpstring = 
	"""

	Usage: julia init [-parameters]

	Parameters:
  	  -f      [data file path]
  	  -N      [number particles]
  	  -i      [initial condition string]


	"""
	print(helpstring)
	return helpstring
end

function unpack_ics(xy,vxy)
	##
	N = length(xy[:,1])
	k = 1:N
	x, y = xy[:,1], xy[:,2]
	vx, vy = vxy[:,1], vxy[:,2]
	##
	return k, x, y, vx, vy
end

function headstring(io,N)
	## print
	for dat in ("NOTES", "N", "r0", "dT", "dt", "dtt", "t")
		print(io, qp(dat))
	end
	for k in 1:N
		for dat in ("k", "x", "y", "vx", "vy")
			print(io, qp(dat))
		end
	end
	print(io,"\n")
end

function datastring(io, N, k, x, y, vx, vy)
	for dat in ("init", N, 0.0, 0.0, 0.0, 0.0, 0.0)
		print(io, qp(dat))
	end
	for k in 1:N
		for dat in (k, x[k], y[k], vx[k], vy[k])
			print(io, qp(dat))
		end
	end
	print(io,"\n")
end

function main(;fname="data.txt", N=3, ic=random)
	## initial data
	data = unpack_ics(ic(N)...)
	## go
	open(fname, "w") do io
		headstring(io, N)
		datastring(io, N, data...)
	end
end

############ run main #############
## defaults if no input
FNAME = "data.txt"
NX = 3
IC = "corner"
## get args
args = ARGS
## check for parameters
for opt in 1:length(args)
	if args[opt] in ("-h",)
		help()
	end
	if args[opt] in ("-f", "-file")
		global FNAME = args[opt+1]
	end
	if args[opt] in ("-N",)
		global NX = parse(Int, args[opt+1])
	end
	if args[opt] in ("-i", "-ic")
		global IC = args[opt+1]
	end
end
## ic string to function
IC = getfield(Main, Symbol(IC))
## run
main(;fname=FNAME, N=1*NX, ic=IC)
###################################

