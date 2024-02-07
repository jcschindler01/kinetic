
# using Revise

# include("../../core/kinetic.jl")

using .Kinetic

dat = Datapoint(N=100, r0=.005, dt=.01, div=10, ic="corner", integrator="free")

T = .05
nsteps = floor(T/dat.dt)
xybins = 3
fbins = 50

function logmult2(N, nvec)
	return logfac2(N) - sum(logfac2.(nvec))
end

function log2_qf(f; N=1000, df=.01)
	r"""
	Compute probability qf to be in macrostate described by coarse-fraction vector f.

	f = 1d vector of fractional values
	b = total number of (equally probable) spatial bins = length(f)
	N = number of particles
	df = coarse fraction increment

	Returns a lower and upper approximation and their midpoint,
	of output value S = log2(qf).

	Theory:
	$q_f = \sum_{\vec{n}} N_{\vec{n}} b^{-N}$.
	"""
	##
	b = length(f)
	##
	eps = 1e-15
	nmin = ceil.(f .* N)					## minimum allowed particles per bin given f
	nmax = floor.((f.+df.-eps) .* N)		## maximum allowed particles per bin given f
	Nmin, Nmax = sum(nmin), sum(nmax)
	##
	Smin, Smax = NaN, NaN
	##
	if (Nmin>N) | (Nmax<N)
		Smin, Smax = -Inf, -Inf
	elseif Nmin==N
		Smax = - N * log2(b) + logmult2(N, nmin)
		Smin = Smax
	elseif Nmax==N
		Smin = - N * log2(b) + logmult2(N, nmax)
		Smax = Smin
	else
		##
		nnmin = 1 * nmin
		nnmax = 1 * nmax
		excess_min = N .- (Nmax .- nmax)	## particles left when all other bins are full
		excess_max = N .- (Nmin .- nmin)	## particles left when all other bins are empty
		nnmin .= max.(nmin,excess_min)
		nnmax .= max.(nmax,excess_max)
		##
		Smin = - N * log2(b) + logmult2(N, nnmax)
		Smax = - N * log2(b) + logmult2(N, nnmin) #+ b * log2(MathConstants.e)
	end
	##
	midpoint = (Smax + Smin)/2
	bounds = (Smin, Smax)
	##
	return midpoint, bounds
end

N = 1000
fbins = 50
f = rand(9)
f = round.(f./sum(f); digits=3)
println(f)

s = log2_qf(f; N=N, df=1/fbins)

println(s)


# vol = qref_spatial(N=dat.N, xybins=xybins, fbins=fbins, trials=1000)

# println(qp(vol))
# println()

# # S0 = Stau(N=dat.N)

# # println(S0)
# # println()

# for k=1:nsteps
# 	evolve!(dat)
# 	ms = macrostate_spatial(dat.xy[:,1], dat.xy[:,2], xybins=xybins, fbins=fbins)
# 	S = S0 + 
# 	println(ms)
# 	println(S0-S)
# end

