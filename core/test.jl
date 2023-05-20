
## load package
include("kinetic.jl")
using .Kinetic

@time begin

# data = Datapoint(N=200, r0=.005, div=10, integrator="naive")
# evolve!(data; tprime=0)

data = Datapoint(N=2000, r0=.005, div=100, integrator="naive")

new_file(data)

t = range(0,.5,61)

for i in 1:length(t)
	evolve!(data; tprime=t[i])
	to_file(data)
	println(qp((i, data.cc)))
end

end

run(`python3 animate.py`)

run(`viewnior out.gif`)
