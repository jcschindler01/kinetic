




include("helpers.jl")
using SpecialFunctions: erf

"""
Sample frequency coarse-grainings

An important class of coarse-grainings are those based on distributions of 
single-particle properties. Suppose H = o H_k and let M be a local POVM in each
subsystem. Then N = o M is a POVM whose outcomes are strings of M outcomes. 
From this one forms the coarser NN = LN, which groups together strings based on
the sample frequency of each outcome. The sample frequencies are taken coarsely
with width df, so that NN_f is the sum over all M strings with sample frequency n/N
in the range [f,f+df). Thus macrostate labels are in general subnormalized, in the 
range 1 - m df <= sum f <= 1.

We can say that tau is iid wrt to M if tr(oM tau) = prod_k tr(M tau_k) = prod_k Q_k
with all tau_k identical. Let q = tr(NN_f tau). If tau is iid with respect to M, then
by Sanov's theorem
-(1/K) log q =~ D(P*||Q)
where Q is the true iid distribution and P* is the worst [largest D(P*||Q)] 
sample distribution among the set allowed by N_f.

We therefore can evaluate macrostate volumes for coarse sample frequency
whenever tau is iid wrt to M, obtaining the value
- log q = K D(P*||Q).
"""

function Stau(;N=1000, alpha=1111, d=2)
	"""
	Entropy of energy-canonical spatial-microcanonical tau.
	Alpha = L/lamtherm.
	"""
	return N * d * log2(alpha*sqrt(ee)) - logfac2(N)
end

function PSTAR(f, df, Q)
	"""
	Start with actual distribution f, not coarse.
	Both f and Q have correct norm and energy automatically.
	Thus their difference is norm and energy preserving.
	Simply apply the biggest multiple of their difference allowed by df/2 range.
	"""
	##
	X = Q - f
	Xscale = maximum(abs.(X))
	dX = min(Xscale, df/2).* (X/Xscale)
	P = f + dX
	##
	return P
end

function log2pp_qf(f, df, M=:spatial; Margs...)
	"""
	Calculate log2(q_f) per particle.

	P is worst case distribution compatible with f,df.

	Q is reference distribution computed from standard tau.
	M is the function to compute Q distribution ie one-particle measurement on tau.
	Margs are parameters of the one-particle measurement.
	"""
	# reference one-particle distribution
	Q = eval(M)(Margs)
	## best actual one-particle distribution
	P = PSTAR(f, df, Q)
	## return relative entropy
	return -D(P, Q)
end

"""
Spatial distribution coarse-graining.
"""

function spatial(args)
	m = args[:xybins]^2
	return ones(m)/m
end

function f_spatial(dat; xybins=6)
	##
	f = hist2d(dat.xy[:,1], dat.xy[:,2], nbins=xybins)/dat.N
	return f[:]
 end

function S_spatial(dat; xybins=6, df=.02)
	S0 = Stau(N=dat.N)
	f = f_spatial(dat; xybins=xybins)
	SM = dat.N * log2pp_qf(f, df, :spatial; xybins=xybins)
	return S0 + SM
end


"""
Speed distribution coarse-graining.
"""

function velocity(args)
	"""
	The single particle speed distribution is obtained from tau by
	changing to spherical velocity coords of each particle, and 
	applying scale conversions from the manuscript. We only need the
	unnormalized behavior (normalize after discretizing) which is
		P(v) ~ v^{d-1} exp(-(1/2)(v/s)^2) dv
	with
		s = sqrt(2E/Nd).
	The indefinite integral for d=2 is proportional to
		-exp(-(1/2)(v/s)^2)
	which means the definite integral from v to v' is
		exp(-(1/2)(v/s)^2) - exp(-(1/2)(v'/s)^2),
	giving the discrete distribution once normalized.		
	"""
	##
	vedges = args[:vedges]
	sigma = args[:sigma]
	v1 = vedges[1:end-1]
	v2 = vedges[2:end]
	vv = (v1+v2)/2
	s = sigma
	Q = exp.(-0.5.*(v1./s).^2) - exp.(-0.5.*(v2./s).^2)
	Q = Q./sum(Q)
	##
	return Q
end

function f_velocity(dat; vedges=0:.1:10)
	##
	f = hist(speeds(dat); nbins=length(vedges)-1, xmin=0, xmax=vedges[end])/dat.N
	return f
 end

## entropy
function S_velocity(dat; dv=.1, min_vmax=5, df=.02)
	##
	S0 = Stau(N=dat.N)
	vmax = max(min_vmax, maximum(speeds(dat))+dv)
	vedges = 0:dv:vmax
	f = f_velocity(dat; vedges=vedges)
	SM = dat.N * log2pp_qf(f, df, :velocity; vedges=vedges, sigma=sigma(dat))
	##
	return S0 + SM
end

"""
Velocity coarse-graining
(vector)
"""

## entropy
function S_vvector(dat; dv=.5, min_vmax=1, df=.02)
	##
	S0 = Stau(N=dat.N)
	nbins = ceil(Int,min_vmax/dv)
	vmax = nbins*dv
	vedges = -vmax:dv:vmax
	f = hist2dv(dat.vxy[:,1], dat.vxy[:,2], xymin=-vmax, xymax=vmax, nbins=2*nbins)/dat.N
	Q = Q_vvector(vedges=vedges, s=sigma(dat))
	Pstar = PSTAR(f, df, Q)
	S = S0 - dat.N * D(Pstar,Q)
	##
	return 1*S
end

function Q_vvector(;vedges=[-3,3,.3], s=1)
	## get prob in each bin for exp(-v^2/2s^2), s=sqrt(2E/Nd)
	integral = erf.(vedges/(sqrt(2)*s))
	q = integral[2:end]-integral[1:end-1]
	qq = q .* q'
	Q = qq/sum(qq)
	return Q
end


"""
Local energy coarse-graining.
"""

function systemA(dat)
	##
	k = 1:dat.N
	mask = k .< dat.N/2
	##
	if dat.ic=="hotcold"
		mask = dat.xy0[:,1] .< 0.5
	end
	if dat.ic=="corner"
		mask = k .<= dat.N/4
	end
	if dat.ic=="gun"
		mask = dat.xy0[:,1] .> 0.7
	end
	if dat.ic=="chain"
		mask = k .< dat.N/4
	end
	return mask
end


## entropy
function S_localE(dat; fA=.5, d=2)
	##
	S0 = Stau(N=dat.N)
	N  = dat.N
	##
	maskA = systemA(dat)
	##
	eAB = dat.vxy[:,1].^2 + dat.vxy[:,2].^2
	eA = eAB[maskA]
	eB = eAB[.!maskA]
	NA = length(eA)
	NB = length(eB)
	##
	E, EA, EB = sum(eAB), sum(eA), sum(eB)
	##
	eeA = NA * E / N
	eeB = NB * E / N
	##
	aA = NA*d/2 - 1
	aB = NB*d/2 - 1
	##
	negent = aA * log2(eeA/EA) + aB * log2(eeB/EB)
	## return
	return S0 - negent
end

function S_localEA(dat; fA=.5, d=2)
	##
	S0 = Stau(N=dat.N)
	N  = dat.N
	NA = floor(Int, fA * N)
	NB = N - NA
	##
	eAB = dat.vxy[:,1].^2 + dat.vxy[:,2].^2
	eA = eAB[1:NA]
	eB = eAB[NA+1:end]
	##
	E, EA, EB = sum(eAB), sum(eA), sum(eB)
	##
	eeA = NA * E / N
	eeB = NB * E / N
	##
	aA = NA*d/2 - 1
	aB = NB*d/2 - 1
	##
	beta = N*d /(2*E)
	##
	negent = aA * log2(eeA/EA) + beta*(EA - eeA)*log2(exp(1))
	##
	## to make it the "conventional OE" version uncomment below
	#### negent = aA * log2(eeA/EA)
	##
	## return
	return S0 - negent
end


