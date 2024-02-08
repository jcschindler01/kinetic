

# include("../../core/kinetic.jl")

using .Kinetic

dat = Datapoint(N=10000, r0=.005, dt=.01, div=10, ic="corner", integrator="free")

T = 1
nsteps = floor(T/dat.dt)
xybins = 5
fbins = 100

S0 = Stau(N=dat.N)

for k=1:nsteps
	evolve!(dat)
	ms = spatial_macrostate(dat.xy[:,1], dat.xy[:,2], xybins=xybins, fbins=fbins, mode="fxy_coarse")
	f = vcat(ms...)
	S = spatial_log2_qf(f; N=dat.N, df=1/fbins)
	println(round.(f; digits=3))
	println(S0/dat.N)
	println(S/dat.N)
	println()
end

