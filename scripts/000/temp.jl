
include("../../core/kinetic.jl")
include("../../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot


## params
T = 4
rate = 1
loops = 1
dat0 = Datapoint(N=1, r0=.002, dt=.01, div=1, ic="random", integrator="free")

## go
bp = Boxplot()
display(bp.fig)
showinit!(bp, dat0)

## number steps
nsteps = floor(T/dat0.dt)

## go
@time for loop=1:loops

	dat = deepcopy(dat0)
	update!(bp,dat)

	for k=1:nsteps
		evolve!(dat)
		update!(bp,dat)
		sleep(dat.dt/rate)
	end

	reversal!(dat)
	update!(bp,dat)

	for k=1:nsteps
		evolve!(dat)
		update!(bp,dat)
		sleep(dat.dt/rate)
	end

	reversal!(dat)
	update!(bp,dat)
end

