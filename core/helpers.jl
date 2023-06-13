
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
maxboltz(v; T=1) = (v./T).*exp.(-v.^2 ./ (2*T))

