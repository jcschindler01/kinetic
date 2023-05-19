
include("quickprint.jl")

length(ARGS)>0 ? fname=ARGS[1] : fname="data.txt"


function freverse(datafile="data.txt")
	## format is ("note", "N", "r0", "it", dT", "dt", "dtt", "t") + N * ("k", "x", "y", "vx", "vy")
	## read
	datastrings = []
	open(datafile, "r") do io
		datastrings = [strip(s) for s in split(chomp(last(readlines(io))), ",", keepempty=true)]
	end
	## process
	note = "frev"
	it = parse(Int, datastrings[8])
	N = Integer(length(datastrings[9:end-1])/5)
	r0 = parse(Float64, datastrings[3])
	dT, dt, dtt, t = (parse(Float64, s) for s in datastrings[4:7])
	x, y, vx, vy = zeros(Float64, N), zeros(Float64, N), zeros(Float64, N), zeros(Float64, N)
	for k in 1:N
		 x[k] = parse(Float64, datastrings[4+5k+1])
		 y[k] = parse(Float64, datastrings[4+5k+2])
		vx[k] = parse(Float64, datastrings[4+5k+3])
		vy[k] = parse(Float64, datastrings[4+5k+4])
	end
	## data
	data = note, N, dT, dt, dtt, t, x, y, vx, vy, it
	## reverse velocities
	vx = -vx
	vy = -vy
	## print
	open(datafile, "a") do io
		for dat in (note, N, r0, dT, dt, dtt, t, it)
			print(io, qp(dat))
		end
		for k in 1:N
			for dat in (k, x[k], y[k], vx[k], vy[k])
				print(io, qp(dat))
			end
		end
		print(io, "\n")
	end
	## return
	return data
end




