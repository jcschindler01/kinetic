
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
phasedist(d1::Datapoint, d2::Datapoint) = norm(hcat(d2.xy.-d1.xy,d2.vxy.-d1.vxy))

## particle pairs
ppairs(N) = [(i,j) for j=1:N, i=1:N if j>i]

## area
area(dat::Datapoint) = dat.N * pi * dat.r0^2

## value parser
parser(T::Type, s::String) = T==String ? s : parse(T,s)

## get initial conditions
getic(ic::String,N::Int) = eval(:($(Symbol(ic))($(N))))


