
include("../../core/kinetic.jl")

using .Kinetic


dat = Datapoint(N=500, r0=.01, dt=.01, div=10, ic="corner", integrator="naive")

S = S_vvector(dat)

println(S)




