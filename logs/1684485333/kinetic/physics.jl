
#= 
Proper way to go from x,y to xy and back.
	xy  = hcat(x,y)
	x,y = xy[:,1], xy[:,2]
The joint one has the form
	xy[N,2]
so that
	xy[k,:]
is the vector position of kth particle.
=#

function timestep(r0, dt, t, x, y, vx, vy, it, pairs)
	# free half step
	x  .+=  vx .* dt/2
	y  .+=  vy .* dt/2
	# update velocities
	velupdate!(r0, dt, t, x, y, vx, vy, it, pairs)
	# free half step
	x  .+=  vx .* dt/2
	y  .+=  vy .* dt/2
	# update time and iterflag
	t  += dt
	it = (it+1)%2
	## return
	return  t, x, y, vx, vy, it
end

function velupdate!(r0, dt, t, x, y, vx, vy, it, pairs)
	## process collisions
	collide!(r0, dt, t, x, y, vx, vy, it, pairs)
end

function reverse!(x, y, vx, vy, it)
	vx = -vx
	vy = -vy
	return  x, y, vx, vy, it
end

function PAIRS(N)
	return [(i,j) for i=1:N, j=1:N if j>i]
end

function collide!(r0, dt, t, x, y, vx, vy, it, pairs; interact=true)
	##
	N = length(x)
	## wall
	wall!(r0, dt, t, x, y, vx, vy, it)
	## interparticle
	if interact==true
		ball!(r0, dt, t, x, y, vx, vy, it, pairs)
	end
end

function wall!(r0, dt, t, x, y, vx, vy, it)
	N = length(x)
	for i=1:N
		## x walls
		if x[i]<0 || x[i]>1
			vx[i] *= -1
		end
		## y walls
		if y[i]<0 || y[i]>1
			vy[i] *= -1
		end
	end
end


function ball!(r0, dt, t, x, y, vx, vy, it, pairs)
	## which way to iterate
	it==0 ? (loop = pairs) : (loop = reverse(pairs))
	## loop over a list of unique (i,j) pairs
	for (i,j) in loop
		## only interact if both in bounds
		if !(0<x[i]<1 && 0<y[i]<1 && 0<x[j]<1 && 0<y[j]<1)
			continue
		end
		## calculate relative position
		dx = 0.5 * (x[j]-x[i])
		dy = 0.5 * (y[j]-y[i])
		dr  = sqrt(dx^2 + dy^2)
		## only interact if sufficiently close
		if !(dr <= 2*r0)
			continue
		end
		## calculate relative velocity
		u  = 0.5*(vx[j]-vx[i])
		v  = 0.5*(vy[j]-vy[i])
		## calculate approach velocity
		eps = 1e-12
		av = (u*dx+v*dy)/(dr+eps)
		## only interact if sufficient magnitude approach velocity
		vmin = 0.0*r0
		if !(abs(av) > vmin)
			continue
		end 
		## calculate center of mass velocity
		U  = 0.5*(vx[j]+vx[i])
		V  = 0.5*(vy[j]+vy[i])
		## calculate new relative velocity
		uprime = u - 2.0*dx*(u*dx+v*dy)/(dx*dx+dy*dy)
		vprime = v - 2.0*dy*(u*dx+v*dy)/(dx*dx+dy*dy)
		## update velocities
		vx[i]  = U - uprime
		vx[j]  = U + uprime
		vy[i]  = V - vprime
		vy[j]  = V + vprime
	end
end

