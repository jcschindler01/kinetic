
## load package
include("kinetic.jl")
using .Kinetic

# @time begin
	
data = Datapoint()

new_file(data)

evolve!(data; integrator="naive", div=3, dt=1)

to_file(data)

# new_file(data)

# t = range(0,1,11)

# for i in 1:length(t)
# 	evolve!(data; tprime=t[i])
# 	to_file(data)
# 	println(qp((i, data.cc)))
# end

# end


# run(`python3 animate.py`)

# run(`viewnior out.gif`)
