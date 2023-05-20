
## random xy and vxy
zrand(N) = round.(rand(Float64, (N,2)); digits=5)
vrand(N) = round.(2. .* rand(Float64, (N,2)) .- 1; digits=5)

## distances
dot(a,b) = a' * b
norm(a) = sqrt(dot(a,a))
dist(a,b) = norm(a-b)
phasedist(d1::Datapoint, d2::Datapoint) = sqrt(dist(d1.xy,d2.xy)^2 + dist(d1.vxy,d2.vxy)^2)
