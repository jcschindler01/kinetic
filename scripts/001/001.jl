
# using Revise

# include("../../core/kinetic.jl")

using .Kinetic

dat = Datapoint(N=2500, r0=.005, dt=.01, div=10, ic="corner", integrator="naive")

T = 4
nsteps = floor(T/dat.dt)
xybins = 2
fbins = 17

vol = qref_spatial(N=dat.N, xybins=xybins, fbins=fbins, trials=100000)

println(qp(vol))
println()

S0 = Stau(N=dat.N)

println(S0)
println()

for k=1:nsteps
	evolve!(dat)
	ms = macrostate_spatial(dat.xy[:,1], dat.xy[:,2], xybins=xybins, fbins=fbins)
	S = S0 + log2(qref(vol,ms))
	println(ms)
	println(S0-S)
end

