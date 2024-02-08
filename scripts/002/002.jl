

# include("../../core/kinetic.jl")

using .Kinetic

dat = Datapoint(N=500, r0=.005, dt=.01, div=10, ic="corner", integrator="free")

T = 1
nsteps = floor(T/dat.dt)
xybins = 3
fbins = 50

function logmult2(N, nvec)
	return logfac2(N) - sum(logfac2.(nvec))
end

function log2_qf(f; N=1000, df=.1)
	"""Approximate method."""
	##
	S = NaN
	b = length(f)
	##
	Df = (1-sum(f))/b					## excess fraction per bin
	dm = round(Int,min(Df,df-Df)*N)		## bin bounds for m value
	nbar = round.((f.+Df)*N)
	##
	if dm<0
		S = -Inf
	elseif dm==0
		S = - N * log2(b) + logmult2(N, nbar)
	elseif dm>0
		S = - N * log2(b) + logmult2(N, nbar) + b * log2(2*dm)
	end
	##
	return S
end



S0 = Stau(N=dat.N)

println(S0)
println()

for k=1:nsteps
	evolve!(dat)
	ms = macrostate_spatial(dat.xy[:,1], dat.xy[:,2], xybins=xybins, fbins=fbins, mode="fxy_coarse")
	f = vcat(ms...)
	S = log2_qf(f; N=dat.N, df=1/fbins)
	println(round.(f; digits=3))
	println(S)
	println()
end


