

include("./ics.jl")
include("./physics.jl")

## params
N = 3
ic = corner
dt = .1
it = 0
r0 = .001
steps = 1

## 
xy, vxy = ic(N)

##
x0,   y0 =  xy[:,1],  xy[:,2]
vx0, vy0 = vxy[:,1], vxy[:,2]

timestep(x, y, vx, vy, it) = timestep(r0, dt, 0, x, y, vx, vy, it, PAIRS(N))

function dance(x, y, vx, vy, it; steps=steps)
	x0, y0, vx0, vy0 = 1*x, 1*y, 1*vx, 1*vy
	for step=1:steps
		t, x, y, vx, vy, it = timestep(x, y, vx, vy, it)
	end
	reverse!(x, y, vx, vy, it)
	for step=1:steps
		t, x, y, vx, vy, it = timestep(x, y, vx, vy, it)
	end
	println()
	println(sum(abs.(x-x0)))
	println(sum(abs.(y-y0)))
	println(sum(abs.(vx-vx0)))
	println(sum(abs.(vy-vy0)))
	println()
end



dance(x0, y0, vx0, vy0, it)
