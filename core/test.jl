
## load package
include("kinetic.jl")
using .Kinetic

data = Datapoint(N=100, r0=.005, div=1000, integrator="naive")

#new_file(data)


tt = range(0,4,101)

@time for i in 1:length(tt)
	println(qp((i, data.cc)))
	evolve!(data; tprime=tt[i])
	#to_file(data)
end



# run(`python3 animate.py`)

# run(`viewnior out.gif`)
