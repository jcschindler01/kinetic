

include("./ics.jl")
include("./physics.jl")
include("./quickprint.jl")

## params
N = 3
ic = circle
dt = .001
it = 0
r0 = .01
steps = 4000

## 
xy, vxy = ic(N)

##
x0,   y0 =  xy[:,1],  xy[:,2]
vx0, vy0 = vxy[:,1], vxy[:,2]

timestep(x, y, vx, vy, it) = timestep(r0, dt, 0, x, y, vx, vy, it, PAIRS(N))[2:end]

trev(x, y, vx, vy, it) = x, y, -vx, -vy, it

function pp(data; p=2)
	for k in p:p
		print(qp(k))
		for dat in data[1:end-1]
			print(qp(dat[k]))
		end
		print(qp(data[end]))
	end
	print("\n")
end

function dance(x, y, vx, vy, it; steps=steps)
	## dance
	data0 = x, y, vx, vy, it
	data  = 1 .* data0
	pp(data0)
	for i in 1:steps
		data = timestep(data...)
		pp(data)
	end
	data = trev(data...)
	pp(data)
	for i in 1:steps
		data = timestep(data...)
		pp(data)
	end
	## eval
	diff = data .- trev(data0...)
	q = sum(collect(sum(abs.(diff[k])) for k in 1:length(data)))
	println()
	println("q = ", q)
	println()
	for k=1:N
		pp(diff; p=k)
	end
	println()
	println("q = ", q)
	println()
end



dance(x0, y0, vx0, vy0, it)
