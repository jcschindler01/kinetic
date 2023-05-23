
## load package
include("kinetic.jl")
using .Kinetic


@time begin

data = Datapoint(N=500, r0=.005, div=100, integrator="sym")

evolve!(data; tprime=4, div=1000)
	
new_file(data)

t = range(0,.5,11)

for i in 1:length(t)
	evolve!(data; tprime=t[i], div=100)
	to_file(data)
	println(qp((i, data.cc)))
end

end


run(`python3 animate.py`)

run(`viewnior out.gif`)
