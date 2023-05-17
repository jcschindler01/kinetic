
# Usage: julia kinetic.jl [file] [-parameters]

include("./quickprint.jl")
include("./physics.jl")


function help()
	helpstring = 
	"""

	Usage: julia kinetic.jl [file] [-parameters]

	Parameters:
	  -dT     [propagation time]
	  -dt     [output timestep]
	  -dtt    [internal timestep]
	  -nsteps [total output timesteps] (nsteps,dT overrides dt)
  	  -div    [internal timesteps per dt] (div,dt overrides dtt)

	"""
	print(helpstring)
	return helpstring
end

function integify(x; eps=1e-9)
	if abs(x - round(x))< eps
		return Integer(round(x))
	else
		return nothing
	end	
end

function is_zeroish(x; precis=6)
	return round(x; digits=precis)==0
end

struct Params
	file::String         # data file
	dT::Real             # propagation time
	dt::Real             # output timestep
	dtt::Real            # internal timestep
	nsteps::Integer	     # total output timesteps dT/dt
	div::Integer         # internal timesteps per output dt/dtt
	NSTEPS::Integer      # total internal timesteps dT/dtt
	r0::Real             # hard sphere radius
end


function Params(;file="data.txt", dT=1.0, dt=1e-1, dtt=1e-1, div=nothing, nsteps=nothing, r0=.01)
	# Uses default values unless given in args. dt = dT/nsteps. dtt=dt/div.
	if nsteps != nothing
		dt = dT/nsteps
	end
	if div != nothing
		dtt = dt/div
	end
	nsteps = integify(dT/dt)
	div = integify(dt/dtt)
	NSTEPS = integify(dT/dtt)
	return Params(file, dT, dt, dtt, nsteps, div, NSTEPS, r0)
end


function readdata(datafile)
	## format is ("note", "N", "r0", "dT", "dt", "dtt", "t") + N * ("k", "x", "y", "vx", "vy")
	## read
	datastrings = []
	open(datafile, "r") do io
		datastrings = [strip(s) for s in split(chomp(last(readlines(io))), ",", keepempty=true)]
	end
	## process
	note = datastrings[1]
	N = Integer(length(datastrings[8:end-1])/5)
	dT, dt, dtt, t = (parse(Float64, s) for s in datastrings[4:7])
	x, y, vx, vy = zeros(Float64, N), zeros(Float64, N), zeros(Float64, N), zeros(Float64, N)
	for k in 1:N
		 x[k] = parse(Float64, datastrings[3+5k+1])
		 y[k] = parse(Float64, datastrings[3+5k+2])
		vx[k] = parse(Float64, datastrings[3+5k+3])
		vy[k] = parse(Float64, datastrings[3+5k+4])
	end
	## data
	data = note, N, dT, dt, dtt, t, x, y, vx, vy
	## return
	return data
end

function simulate(params, initial_data)
	## unpack initial data
	note0, N0, dT0, dt0, dtt0, t0, x0, y0, vx0, vy0 = initial_data
	N, dT, dt, dtt, r0 = N0, params.dT, params.dt, params.dtt, params.r0
	## go
	T = 0
	nstep = 1
	note = "sim"
	t, x, y, vx, vy = t0, x0, y0, vx0, vy0
	open(params.file, "a") do io
		while T < dT
			T += dtt
			t, x, y, vx, vy = timestep(r0, dtt, t, x, y, vx, vy)
			if is_zeroish(T-nstep*dt)
				println(nstep)
				nstep +=1
				## print
				for dat in (note, N, r0, dT, dt, dtt, t)
					print(io, qp(dat))
				end
				for k in 1:N
					for dat in (k, x[k], y[k], vx[k], vy[k])
						print(io, qp(dat))
					end
				end
				print(io, "\n")
				note = "..."
				## end print
			end
		end
	end
end


function timestep(r0, dt, t, x, y, vx, vy)
	# xy(t+dt)
	x  .+=  vx .* dt
	y  .+=  vy .* dt
	# vxy(t+dt)
	vx .+= 0
	vy .+= 0
	# t
	t += dt
	## return
	return  t, x, y, vx, vy
end


####### main ########
function main(params)
	## read data from file
	initial_data = readdata(params.file)
	## do simulation
	@time simulate(params, initial_data)
	## return
	return nothing
end
#####################


############ run main #############
## defaults if no input
FNAME="data.txt"
dT=1.0
dt=nothing
dtt=nothing
DIV=10
NSTEPS=10
R0=.01
## get args
args = ARGS
## get file
if length(args)!=0 && args[1][1]!='-'
	FNAME = args[1]
end
## check for params
for opt in 1:length(args)
	if args[opt] in ("-h",)
		help()
	end
	if args[opt] in ("-T", "-dT")
		global dT = parse(Float64, args[opt+1])
	end
	if args[opt] in ("-dt",)
		global dt = parse(Float64, args[opt+1])
	end
	if args[opt] in ("-dtt",)
		global dtt = parse(Float64, args[opt+1])
	end
	if args[opt] in ("-d", "-div")
		global DIV = parse(Float64, args[opt+1])
	end
	if args[opt] in ("-n", "-nsteps")
		global NSTEPS = parse(Float64, args[opt+1])
	end
	if args[opt] in ("-r", "-r0")
		global R0 = parse(Float64, args[opt+1])
	end
end
## params
params = Params(;file=FNAME, dT=dT, dt=dt, dtt=dtt, div=DIV, nsteps=NSTEPS, r0=R0)
## run
main(params)
###################################




