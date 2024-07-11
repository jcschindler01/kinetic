
include("../../core/kinetic.jl")

using .Kinetic


dat = Datapoint(N=500, r0=.01, dt=.01, div=10, ic="gun", integrator="naive")

S = S_localE(dat)

println(S)




