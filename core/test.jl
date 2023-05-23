
## load package
include("kinetic.jl")
using .Kinetic


data = datapoint(N=5, div=1, integrator="sym")

evolve!(data)


# @time begin

# data = Datapoint(N=5, r0=.01, div=1, integrator="sym")


	
# new_file(data)

# t = range(0,1,21)

# for i in 1:length(t)
# 	evolve!(data; tprime=t[i])
# 	to_file(data)
# 	println(qp((i, data.cc)))
# end

end


# run(`python3 animate.py`)

# run(`viewnior out.gif`)
