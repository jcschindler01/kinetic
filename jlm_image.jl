"""
To generate system image:

>> julia

using PackageCompiler

create_sysimage(
	["GLMakie"]; 
	sysimage_path="jlm.so",
	precompile_execution_file="jlm_image.jl",
	)

"""

include("core/kinetic.jl")
include("core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot


## params
T = 2
rate = 2
loops = 1
dat0 = Datapoint(N=50, r0=.005, dt=.01, div=1, ic="corner", integrator="naive")

## go
bp = Boxplot()
display(bp.fig)

## go
for loop=1:loops

    dat = deepcopy(dat0)
    update!(bp,dat)

    while dat.t<=T
        evolve!(dat)
        update!(bp,dat)
        sleep(dat.dt/rate)
    end

    reversal!(dat)

    while dat.t<=2*T
        evolve!(dat)
        update!(bp,dat)
        sleep(dat.dt/rate)
    end
end


