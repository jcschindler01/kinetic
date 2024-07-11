
#=
Back and forth.
	xy  = hcat(x,y)
	x,y = xy[:,1], xy[:,2].
=#


## random xy and vxy
zrand(N) = round.(rand(Float64, (N,2)); digits=5)
vrand(N) = round.(2. .* rand(Float64, (N,2)) .- 1; digits=5)
rms(vxy) = sqrt(sum(vxy.^2)/length(vxy[:,1]))
vrmsnorm(vxy) = round.(vxy/rms(vxy); digits=5)

## distances
dot(a,b) = a' * b
using LinearAlgebra: norm
phasedist(d1, d2) = norm(hcat(d2.xy.-d1.xy,d2.vxy.-d1.vxy))
phasedist(dat) = norm(hcat(dat.xy.-dat.xy0,dat.vxy.-dat.vxy0))

## particle pairs
ppairs(N) = [(i,j) for j=1:N, i=1:N if j>i]

## area
area(dat) = dat.N * pi * dat.r0^2

## value parser
parser(T::Type, s::String) = T==String ? s : parse(T,s)

## get initial conditions
getic(ic::String,N::Int) = eval(:($(Symbol(ic))($(N))))

## log2 factorial
function logfac2(n)
	out = 0.0
	for k=1:n
		out += log2(k)
	end
	return 1.0*out
end

## log2 multinomial coefficients
logmult2(N, nvec) = logfac2(N) - sum(logfac2.(nvec))

## relative entropy
function D(p,q)
	return sum(plog2q(p,p./q))
end

function plog2q(p,q)
	s = NaN*p
	mask = (p.>0).&&(q.>0)
	s[mask] .= p[mask].*log2.(q[mask])
	mask = (p.==0).&&(isfinite.(q))
	s[mask] .= 0
	return s
end

## using DataStructures

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

## create 2d histogram
function hist2dv(vx, vy; nbins=3, xymin=0, xymax=1, eps=1e-12)
	h = zeros(Int, (nbins,nbins))
	dxy = (xymax-xymin)/nbins
	dd = (xymax-xymin)
	x = vx .- xymin
	y = vy .- xymin
	for i=1:length(x)
		xx = clamp(x[i], eps, dd-eps)
		yy = clamp(y[i], eps, dd-eps)
		h[1+floor(Int,xx/dxy), 1+floor(Int,yy/dxy)]+=1
	end
	return h
end


#=
For 2d ideal gas of N particles.

Z0 = 2piV/Bm  						## 1-particle partition function
Z = Z0^N 							## partition function
U = -d/dB log Z = N/B 				## total energy
p(v) dv = Bm v exp(-Bmv^2/2) dv 	## 1-particle speed distribution

Here m=1 and V=1. Define T=1/B (ie kb=1). Then T=U/N.
=#

## thermo
speeds(dat) = sqrt.(dat.vxy[:,1].^2 + dat.vxy[:,2].^2)
energy(dat) = sum(dat.vxy[:,1].^2 + dat.vxy[:,2].^2)/2
temp(dat) = energy(dat)/dat.N
maxboltz(v; T::Real=1) = (v./T).*exp.(-v.^2 ./ (2*T))
beta(dat) = dat.N/energy(dat)
sigma(dat) = sqrt(energy(dat)/dat.N)

## info
const l2e = log2(MathConstants.e)
const ee = MathConstants.e

## dotsize (r,rpx)=[(.004,5) ... (.02,21)]
dotsize(r) = r>0 ? 1 + 1000*r : 5
