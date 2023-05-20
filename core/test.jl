
## load package
include("kinetic.jl")
using .Kinetic

@time begin

data = Datapoint(N=200, r0=.005, div=10, integrator="naive")
new_file(data)
evolve!(data; tprime=0, div=34, integrator="free")
to_file(data)

# data = Datapoint(N=500, r0=.005, div=100, integrator="naive")

# new_file(data)

# t = range(0,.5,61)

# for i in 1:length(t)
# 	evolve!(data; tprime=t[i])
# 	to_file(data)
# 	println(qp((i, data.cc)))
# end

end

# run(`python3 animate.py`)

# run(`viewnior out.gif`)
