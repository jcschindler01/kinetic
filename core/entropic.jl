




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
function hist(x; nbins=10, xmin=0, xmax=1, eps=1e-12)
	h = zeros(Int, nbins)
	dx = (xmax-xmin)/nbins
	for i=1:length(x)
		xx = clamp(x[i], xmin+eps, xmax-eps)
		h[1+floor(Int,xx/dx)]+=1
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
	"""Entropy of energy-canonical spatially-microcanonical reference state tau.
	The parameter alpha relates to the thermal de Broglie wavelength, derives
	from scales chosen in manuscript."""
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

##
function S_spatial(dat; xybins=3, fbins=50)
	ms = spatial_macrostate(dat.xy[:,1], dat.xy[:,2], xybins=xybins, fbins=fbins, mode="fxy_coarse")
	f = vcat(ms...)
	S0 = Stau(N=dat.N)
	Sm = spatial_log2_qf(f; N=dat.N, df=1/fbins)
	return S0 + Sm
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
function velocity_log2_qf(f; N=1000, df=.1, vedges=0:.1:3, sigma=1)
	"""
	Speed histogram macrostates.

	From unitful params in manuscript:

	beta p^2/2m = (beta m L^2/T^2) v^2/2 = (Nd/2E) v^2/2 (rhs unitless values).

	Let sigma = sqrt(2E/Nd). Then

	dp2/Z = dv2/sqrt(2 pi sigma^2) and beta p^2/2m = v^2/(2 sigma^2).

	Thus

	qf = q_f = tr(M_f tau)
	   = prod_N int dp2/Z^2 exp(-beta p^2/2m) M_f
	   = prod_N int dv2/sqrt(2pi sigma^2) exp(-v^2/2sigma^2) M_f
	   = prod_N int sqrt(2pi) v dv/sigma exp(-v^2/2sigma^2) M_f.

	Bin the velocities in phase space. Each point in phase space gives a 
	string of velocity bins. Each string with same n has the same probability.
	Phase space integral becomes sum over strings. M_f projects onto strings
	with compatible n for f. Then

	qf = sum_n N_n prod_N sqrt(2pi) v dv/sigma exp(-v^2/2sigma^2) M_f.

	Again approximate sum_n N_n = (2dm)^b N_nbar. Then

	qf = (2dm)^b N_nbar prod_vbar [sqrt(2pi) v dv/sigma exp(-v^2/2sigma^2)]^nbar

	"""
	##
	S = NaN
	b = length(f)
	dv = vedges[2] - vedges[1]
	vbar = (vedges[1:end-1]+vedges[2:end])/2
	##
	Df = (1-sum(f))/b					## excess fraction per bin
	dm = round(Int,min(Df,df-Df)*N)		## bin bounds for m value
	nbar = round.((f.+Df)*N)
	##
	S_sqbin = b*log2(2*dm)
	S_Nnbar = logmult2(N, nbar)
	S_dv = N*log2(sqrt(2*pi)*dv/sigma)
	SS_prob = nbar .* log2.(vbar .* exp.(-vbar.^2 ./ (2 .* sigma)))
	S_prob = sum(SS_prob)
	total = S_sqbin + S_Nnbar + S_dv + S_prob
	##
	println(S_sqbin)
	println(S_Nnbar)
	println(S_dv)
	println(S_prob)
	println(total)
	println(SS_prob)
	if dm<0
		S = -Inf
	elseif dm==0
		S = S_Nnbar
	elseif dm>0
		S = S_sqbin + S_Nnbar
	end
	##
	return S
end



##
function velocity_macrostate(v; vedges=0:.1:3, fbins=10)
	##
	N = length(v)
	df = 1/fbins
	dv = vedges[2] - vedges[1]
	##
	Nv = hist(v, nbins=length(vedges)-1, xmin=vedges[1], xmax=vedges[end])
	fv = Nv/N
	fv_int = floor.(Int,fv./df)
	fv_coarse = fv_int.*df
	return fv_coarse
end



##
function S_velocity(dat; vedges=0:.1:3, fbins=50)
	f = velocity_macrostate(speeds(dat); vedges=vedges, fbins=fbins)
	S0 = Stau(N=dat.N)
	Sm = velocity_log2_qf(f, N=dat.N, df=1/fbins, vedges=vedges, sigma=sigma(dat))
	return S0 + Sm
end

