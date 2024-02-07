




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

## return qref
function qref(v::Volume, macrostring::String)
	if macrostring in keys(v.vdict)
		return v.vdict[macrostring]/v.Nmc
	else
		return 1/v.Nmc
	end
end

## xy data to spatial density macrostate
function macrostate_spatial(x, y; xybins=3, fbins=10, mode="fxy_string")
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


## create spatial qref
function qref_spatial(;N=100, xybins=3, fbins=10, trials=1000)
	##
	v = Volume()
	v.cg = "spatial density: xybins=$(xybins), fbins=$(fbins)"
	v.tau = "N=$(N), spatially uniform"
	##
	for mc in 1:trials
		##
		x, y = rand(N), rand(N)
		fxy_string = macrostate_spatial(x, y, xybins=xybins, fbins=fbins)
		v.Nmc += 1
		if fxy_string in keys(v.vdict)
			v.vdict[fxy_string] += 1
		else
			v.vdict[fxy_string]  = 1
		end
	end
	sort!(v.vdict, byvalue=true, rev=true)
	##
	return v
end

## quick print
function qp(v::Volume)
	println()
	println(v.cg)
	println(v.tau)
	println(v.Nmc)
	show = collect(keys(v.vdict))
	stop = min(length(show), 10)
	for x in show[1:stop]
		println("$(x) q=$(qref(v,x))")
	end
	if stop!=length(show)
		println("...")
	end
	println()
end

##
function logfac2(n)
	out = 0.0
	for k=1:n
		out += log2(k)
	end
	return 1.0*out
end

##
function Stau(;d=2,N=1000,alpha=1111)
	ee = MathConstants.e
	S = N*d*log2(alpha*sqrt(ee)) - logfac2(N)
	return 1.0*S
end

