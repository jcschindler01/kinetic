

cd(@__DIR__)

include("../core/kinetic.jl")
include("../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot

## params
T = 5
rate = 1
dat = Datapoint(N=500, r0=.005, dt=.005, div=10, ic="hotcold", integrator="naive")

##
save = true

## 
tempfile = "temp.txt"
outfile = "txt/$(dat.ic)_$(dat.integrator)_N$(dat.N)_T$(T)_$(round(Int,time())).txt"

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
if save==true
    cp(tempfile, outfile)
end



