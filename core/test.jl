
## load package
include("kinetic.jl")
using .Kinetic


@time begin

data0 = Datapoint(N=10, r0=.002, div=10, ic="corner", integrator="sym")
data0T = deepcopy(data0)
reversal!(data0T)
data  = deepcopy(data0)
new_file(data)

T = 1
frames = 31

dt = T / frames

## forward
for i in 1:frames
	evolve!(data; dt=dt)
	to_file(data)
	println(qp((i, data.cc, phasedist(data,data0)/data.N)))
end
## flip
println("reversal!")
reversal!(data)
to_file(data)
## backward
for i in 1:frames
	evolve!(data; dt=dt)
	to_file(data)
	println(qp((i, data.cc, phasedist(data,data0T)/data.N)))
end
#

end


# run(`python3 animate.py`)

# run(`viewnior out.gif`)
