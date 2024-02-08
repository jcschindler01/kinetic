

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
