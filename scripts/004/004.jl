
include("../../core/kinetic.jl")

using .Kinetic


## params
dat = Datapoint(N=10000, r0=.005, dt=.01, div=10, ic="random", integrator="free")


vedges = 0:.1:5
fbins = 50

S = S_velocity(dat, vedges=vedges, fbins=fbins)

