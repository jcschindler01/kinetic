
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



##
function gun(N; l=.4, r=.05, bullets=10)
	## target
	sq = ceil(Int, sqrt(N-bullets))
	x0 = l*vcat([i/(sq+1) for i=1:sq, j=1:sq]...)
	y0 = l*vcat([j/(sq+1) for i=1:sq, j=1:sq]...)
	xy0  = hcat(x0,y0)[1:N-bullets,:]
	vxy0 = 0 .*xy0[1:N-bullets,:]
	## bullets
	s = range(1,1.5,bullets)
	x1 = (1-1.5*r) .+ r .* cos.(pi*s)
	y1 = (1-1.5*r) .+ r .* sin.(pi*s)
	xy1  = hcat(x1,y1)
	vxy1 = -1 .+ 0 .* xy1
	## 
	xy  = vcat(xy0,xy1)
	vxy = vcat(vxy0,vxy1)
	##
	return xy, vrmsnorm(vxy)
end


## hot left cold right
function hotcold(N; l=.5, TL=10, TR=1)
	## number
	NL = floor(Int, N/2)
	NR = N - NL
	## left
	xyL  = zrand(NL)
	xyL[:,1] .= l .* xyL[:,1]
	vxyL = TL * randn((NL,2))
	## left
	xyR  = zrand(NR)
	xyR[:,1] .= l .+ (1-l).*xyR[:,1]
	vxyR = TR * randn((NR,2))
	## 
	xy  = vcat(xyL,xyR)
	vxy = vcat(vxyL,vxyR)
	##
	return xy, vrmsnorm(vxy)
end
