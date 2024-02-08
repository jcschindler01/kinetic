




using DataStructures

## volume data type
@kwdef mutable struct Volume
	cg::String = "none"								## coarse graining info
	tau::String = "none"							## reference state info
	Stau::Float64 = 0.								## S(tau)
	Nmc::Int = 0									## number of monte carlo trials
	vdict::OrderedDict{String,Int}= OrderedDict()	## volume lookup dict
end

## create histogram with even bin widths
function hist(x; nbins=10, xmin=0, xmax=1)
	h = zeros(Int, nbins)
	dx = (xmax-xmin)/nbins
	for i=1:length(x)
		h[1+floor(Int,x[i]/dx)]+=1
	end
	return h
end

## create 2d histogram
function hist2d(x, y; nbins=3, xymin=0, xymax=1, eps=1e-12)
	h = zeros(Int, (nbins,nbins))
	dxy = (xymax-xymin)/nbins
	for i=1:length(x)
		xx = clamp(x[i], xymin+eps, xymax-eps)
		yy = clamp(y[i], xymin+eps, xymax-eps)
		h[1+floor(Int,xx/dxy), 1+floor(Int,yy/dxy)]+=1
	end
	return h
end

## xy data to spatial density macrostate
function spatial_macrostate(x, y; xybins=3, fbins=10, mode="fxy_coarse")
	N = length(x)
	df = 1/fbins
	Nxy = hist2d(x, y, nbins=xybins)
	fxy = Nxy/N
	fxy_int = floor.(Int,fxy./df)
	fxy_coarse = fxy_int.*df
	fxy_string = string(fxy_int)
	out = Dict([("Nxy",Nxy),("fxy",fxy),("fxy_int",fxy_int),("fxy_string", fxy_string),("fxy_coarse",fxy_coarse)])
	out = (mode=="dict" ? out : out[mode])
	return out
end

##
function Stau(;d=2,N=1000,alpha=1111)
	ee = MathConstants.e
	S = N*d*log2(alpha*sqrt(ee)) - logfac2(N)
	return 1.0*S
end

##
function spatial_log2_qf(f; N=1000, df=.1)
	"""
	Spatial density macrostates.

	Approximately calculate log2(q_f) for the macrostate 
	described by coarse fraction f with bin width df,
	with N particles and b spatial bins.

	The exact value is 
	q_f = sum_n N_n b^{-N}

	where the sum runs over all n = (n_1, n_2, ..., n_b) compatible 
	with total particles sum(n)=N and coarse fraction n/N ~ f.
	Here N_n = multinomial(N, n) = N! / prod(n!).

	A coarse fraction f must sum between 1-b*dF and 1 in order to be valid.
	Generally it sums to less than 1. In this case define
	Df = (1 - sum(f))/b
	as the extra fraction per bin.

	So long as 0<=Df<=df, the coarse fraction f and width df are valid and compatible.
	If this is false, no configuration works, so S=-Inf.

	When 0<=Df<=df is true we can explicitly write a working configuration.
	Take nbar = (f+Df)*N, which sums to N (mod rounding) and is in the correct bin.
	We then write the exact value as
	q_f = sum_m N! / prod((nbar+m)!) b^{-N}
	where the sum is over all m = (m_1, ..., m_b) compatible
	with total particles sum(m)=0 and coarse fraction (nbar+m)/N ~ f.

	Now we start approximating in the case where multiple terms contribute.

	Define dm = N*min(Df, df-Df), which in the valid case is between 0 and Df/2,
	and which is the biggest symmetrical range of m values contained in each bin.
	This gives an estimate of how many terms contribute to the sum, because
	when nbar is close to bin edges, not many terms can contribute. But within this
	range, for small numbers of bins, a good fraction of these terms contribute.

	We approximate the sum as having (2*dm)^b equal terms, each evaluated at nbar.

	I don't know exactly how bad the approximation is. 
	It at least has some mutually cancelling factors, in that
	we discount some terms from the full m sum, but count extra in the dm sum.
	"""
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


