
using Printf

## load package
include("kinetic.jl")
using .Kinetic


@time begin

data0 = Datapoint(N=10, r0=.005, div=100, ic="circle", integrator="naive")
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
	println(qp((i, data.cc, data.mcc, phasedist(data,data0)/data.N)))
end
## flip
println("reversal!")
reversal!(data)
to_file(data)
## backward
for i in 1:frames
	evolve!(data; dt=dt)
	to_file(data)
	println(qp((i, data.cc, data.mcc, phasedist(data,data0T)/data.N)))
end
#



println("\nscore =\n$(phasedist(data,data0T))\n")

end


run(`python3 animate.py`)

run(`viewnior out.gif`)
