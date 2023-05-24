
## random
function random(N)
	xy  = zrand(N)
	vxy = vrand(N)
	return xy, vrmsnorm(vxy)
end

## random in corner
function corner(N; l=.1)
	xy  = l .* zrand(N)
	vxy = vrand(N)
	return xy, vrmsnorm(vxy)
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
	return xy, vrmsnorm(vxy)
end

## famous chain number 2
function chain2(N)
	##
	lam = range(0,1,N)
	##
	x = .5 .* lam
	y = x .^ 2
	vx = 1 .+ 0*y
	vy = .4 .+ 0*y
	##
	xy  = hcat(x,y)
	vxy = hcat(vx,vy)
	##
	return xy, vrmsnorm(vxy)
end

## circle
function circle(N; v0=1, R=.25)
	s   = range(0,1,N+1)[1:end-1]
	xy  = hcat(0.5 .+ R*cos.(2*pi*s),0.5 .+ R*sin.(2*pi*s))
	vxy = -v0 .* (xy .- .5)/sum(abs.(xy .- .5))
	return xy, vrmsnorm(vxy)
end

