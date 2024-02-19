




include("helpers.jl")

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

We can approximate P* by distributing the remaining probability optimally.
Let delta = (1 - sum f)/df. Calculate dD = D(P+dP||Q)-D(P||Q). Sort and fill
the first floor(delta) bins with extra df, and bin floor(delta)+1 with extra (delta-mdf)
of probability. This gives an approximate P*. The approximation is that the
df increment could possibly change the dD ordering, which in that case shouldn't
matter too much.

We therefore can evaluate macrostate volumes for coarse sample frequency
whenever tau is iid wrt to M, obtaining the value
- log q = K D(P*||Q).
"""

function Stau(;N=1000, alpha=1111, d=2)
	"""
	Entropy of energy-canonical spatial-microcanonical tau.
	Alpha relates to thermal de Broglie wavelength cf manuscript.
	"""
	return N * d * log2(alpha*sqrt(ee)) - logfac2(N)
end

function PSTAR(f, df, Q)
	##
	m = length(f)
	delta = (1-sum(f))/m
	mfill = floor(Int, delta)
	extra = m*delta - mfill*df
	P = f
	##
	if (delta<0).||(delta>df)
		P = Inf*ones(m)
	else
		## sort most advantageous bins
		dD = plog2q(f.+df, (f.+df)./Q) - plog2q(f, f./Q)
		dD[Q.==0] .= -Inf
		I = sortperm(dD, rev=true)
		## fill most advantageous bins
		P[I[1:mfill]] .+= df
		P[I[mfill+1]] += extra
	end
	##
	return P
end


using Optim

function PSTAR_ENERGY(f, df, Q; e=0.5, dv=.1)
	P0 = f/sum(f)
	m = length(f)
	en = [0.5 .*((n-1)*dv)^2 for n=1:m]
	Df = 1 - sum(f)
	De = e - sum(f.*en)
	emid = De/Df
	fillbins = Df/df
	nlast = findfirst(en .>= 2*emid)
	nfill = floor(Int, fillbins/2)
	dF = zeros(m)
	dF[1:nfill] .+= df
	dF[nlast-nfill:nlast] .+= df
	extra = Df - sum(dF)
	dF[nfill+1] += extra
	P = f .+ dF
	clamp!(P, 0, 1)
	return P0
end

function log2pp_qf(f=ones(3)/4, df=.01, M=:spatial; Margs...)
	"""
	Calculate log2(q_f) per particle.

	P is worst case distribution compatible with f,df.

	Q is reference distribution computed from standard tau.
	M is the function to compute Q distribution ie one-particle measurement on tau.
	Margs are parameters of the one-particle measurement.
	"""
	# reference one-particle distribution
	Q = eval(M)(Margs)
	## actual one-particle distribution
	if M==:spatial
		P = PSTAR(f, df, Q)
	elseif M==:velocity
		e = Dict(Margs)[:sigma]^2
		ve = Dict(Margs)[:vedges]
		dv = ve[2]-ve[1]
		P = PSTAR_ENERGY(f, df, Q; e=e, dv=dv)
	end
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

function f_spatial(dat; xybins=3, fbins=10)
	##
	df = 1/fbins
	f = hist2d(dat.xy[:,1], dat.xy[:,2], nbins=xybins)/dat.N
	fc = div.(f, df) .* df
	return fc[:], df
 end

function S_spatial(dat; xybins=3, fbins=100)
	f, df = f_spatial(dat; xybins=xybins, fbins=fbins)
	S0 = Stau(N=dat.N)
	SM = dat.N * log2pp_qf(f, df, :spatial; xybins=xybins)
	return S0 + SM
end


"""
Speed distribution coarse-graining.
"""

## reference distribution from tau
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
	vedges = args[:vedges]
	sigma = args[:sigma]
	v1 = vedges[1:end-1]
	v2 = vedges[2:end]
	vv = (v1+v2)/2
	s = sigma
	Q = exp.(-0.5.*(v1./s).^2) - exp.(-0.5.*(v2./s).^2)
	Q = Q./sum(Q)
	return Q
end

## coarse sample fraction for data
function f_velocity(dat; vedges=0:.1:10, fbins=10)
	##
	df = 1/fbins
	v = speeds(dat)
	nbins = length(vedges)-1
	f = hist(v; nbins=nbins, xmin=0, xmax=vedges[end])/dat.N
	fc = div.(f, df) .* df
	return fc, df
 end

## entropy
function S_velocity(dat; dv=.1, min_vmax=5, fbins=100)
	vmax = maximum(speeds(dat))+dv
	vmax = max(vmax, min_vmax)
	vedges = 0:dv:vmax
	f, df = f_velocity(dat; vedges=vedges, fbins=fbins)
	S0 = Stau(N=dat.N)
	SM = dat.N * log2pp_qf(f, df, :velocity; vedges=vedges, sigma=sigma(dat))
	return S0 + SM
end




