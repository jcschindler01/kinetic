
## helpers 

zzero(N) = zeros((N,2))
zrand(N) = rand(Float64, (N,2))
vrand(N) = 2. .* rand(Float64, (N,2)) .- 1

function vnorm(vxy)
	N = length(vxy[:,1])
	vrms = sqrt(sum(vxy .^ 2)/N)
	if vrms>0
		return vxy ./ vrms
	else
		return vxy
	end
end

## random
function random(N)
	xy  = zrand(N)
	vxy = vrand(N)
	return xy, vnorm(vxy)
end

## random in corner
function corner(N; l=.05)
	xy  = l .* zrand(N)
	vxy = vrand(N)
	return xy, vnorm(vxy)
end

## famous chain
function chain(N)
	##
	lam = range(0,1,N)
	##
	x = .5 .* lam
	y = x .^ 2
	vx = exp.(3 .* x)
	vy = exp.(3 .* y)
	##
	xy  = hcat(x,y)
	vxy = hcat(vx,vy)
	##
	return xy, vnorm(vxy)
end

## circle
function circle(N; v0=0.00, R=.1)
	s   = range(0,1,N+1)[1:end-1]
	xy  = hcat(0.5 .+ R*cos.(2*pi*s),0.5 .+ R*sin.(2*pi*s))
	vxy = -v0 .* (xy .- .5)
	return xy, vxy
end


