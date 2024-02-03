
include("../../core/kinetic.jl")
include("../../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot


## params
T = 2
rate = 10
loops = 1
dat0 = Datapoint(N=200, r0=.005, dt=.01, div=10, ic="random", integrator="naive")

## go
bp = Boxplot()
display(bp.fig)
showinit!(bp, dat0)

## file
new_file(dat0)

## number steps
nsteps = floor(T/dat0.dt)

## go
for loop=1:loops

	dat = deepcopy(dat0)
	update!(bp,dat)

	for k=1:nsteps
		evolve!(dat)
		to_file(dat)
		update!(bp,dat)
		sleep(dat.dt/rate)
	end

	reversal!(dat)
	update!(bp,dat)

	for k=1:nsteps
		evolve!(dat)
		to_file(dat)
		update!(bp,dat)
		sleep(dat.dt/rate)
	end

	reversal!(dat)
	update!(bp,dat)
end

