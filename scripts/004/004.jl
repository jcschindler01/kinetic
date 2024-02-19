
include("../../core/kinetic.jl")

using .Kinetic


## params
dat = Datapoint(N=500, r0=.005, dt=.01, div=10, ic="gun", integrator="free")


S = S_velocity(dat)

println(S)
