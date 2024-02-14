

cd(@__DIR__)

include("../core/kinetic.jl")
include("../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot


## params
T = .03
rate = 1
dat = Datapoint(N=3, r0=.005, dt=.01, div=10, ic="random", integrator="free")

## 
tempfile = "./_temp.txt"
outfile = "$(dat.ic)_$(dat.integrator)_N$(dat.N)_T$(T)_$(round(Int,time())).txt"

##
new_file(dat, tempfile)
to_file(dat, tempfile)

## number steps
nsteps = floor(T/dat.dt)

## make plot
bp = Boxplot()
display(bp.fig)
showinit!(bp,dat)
update!(bp,dat)

## go
for k=1:nsteps
    evolve!(dat)
    update!(bp,dat)
    to_file(dat, tempfile)
    sleep(dat.dt/rate)
end

## copy file
cp(tempfile, outfile)