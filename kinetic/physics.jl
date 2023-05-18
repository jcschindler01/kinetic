
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


function collide!(r0, dt, t, x, y, vx, vy)
	##
	N = length(x)
	## wall
	for i=1:N
		if x[i]<0 || x[i]>1
			vx[i] *= -1
		end
		if y[i]<0 || y[i]>1
			vy[i] *= -1
		end
	end
	## interparticle
	##
	return nothing
end

