
include("../../core/kinetic.jl")
include("../../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot


## params
T = 4
rate = 1
dat = Datapoint(N=200, r0=.01, dt=.01, div=10, ic="corner", integrator="naive")

## number steps
nsteps = floor(T/dat.dt)

## make plot
bp = Boxplot()
display(bp.fig)
update!(bp,dat)

## go
for k=1:nsteps
	evolve!(dat)
	update!(bp,dat)
	sleep(dat.dt/rate)
end




