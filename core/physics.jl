

function free!(dat::Datapoint, dt)
	## free evolution for time dt
	dat.t    += dt
	dat.xy  .+= dt .* dat.vxy
end

function collide!(dat::Datapoint, i::Int, j::Int)
	#= Reverse the component of relative velocity projected onto relative position. =#
	## relative
	dz = 0.5*( dat.xy[j,:] .-  dat.xy[i,:])
	dv = 0.5*(dat.vxy[j,:] .- dat.vxy[i,:])
	## update
	dat.vxy[i,:] .+= 2.0 * (dot(dv,dz)/dot(dz,dz)) .* dz
	dat.vxy[j,:] .-= 2.0 * (dot(dv,dz)/dot(dz,dz)) .* dz
	## collision counter
	dat.cc += 1
	##
end

function wallv!(dat::Datapoint)
	#= reverse velocity if out of bounds =#
	mask = (dat.xy .<= 0) .|| (dat.xy .>= 1)
	dat.vxy[mask] .*= -1
end

function wallz!(dat::Datapoint)
	#= reverse velocity and reflect position if out of bounds =#
	## lower boundary 0
	mask = (dat.xy .<= 0)
	dat.vxy[mask] .*= -1
	dat.xy[mask]    = -dat.xy[mask]
	## upper boundary 1
	mask = (dat.xy .>= 1)
	dat.vxy[mask] .*= -1
	dat.xy[mask]  .*=  2 .- dat.xy[mask]
end

function reversal!(dat::Datapoint)
	## reverse velocities
	dat.vxy .*= -1
	dat.origin = "reverse"
end

function naivecol!(dat::Datapoint)
	#= Collide if close and approaching. =#
	## faster arrays
	x,  y, vx, vy, r0  =  dat.xy[:,1],  dat.xy[:,2], dat.vxy[:,1], dat.vxy[:,2], 1*dat.r0
	## mulitiparticle collision events
	hasinteracted = Set{Int}([])
	## loop over particle pairs
	for i=1:dat.N
		for j=i+1:dat.N
			## values
			dx  =  x[j] -  x[i]
			dy  =  y[j] -  y[i]
			dvx = vx[j] - vx[i]
			dvy = vy[j] - vy[i]
			## conditions
			isclose = sqrt(dx^2 + dy^2) < 2*r0
			isapproaching = dx*dvx+dy*dvy<0
			inbounds = 0<x[i]<1 && 0<y[i]<1 && 0<x[j]<1 && 0<y[j]<1
			## pass if sufficiently close
			if isclose && isapproaching && inbounds
				collide!(dat, i, j)
				i in hasinteracted || j in hasinteracted ? dat.mcc+=1 : nothing
				push!(hasinteracted, i, j)
			end
		end
	end
end

function symcol!(dat::Datapoint)
	#= Time symmetric collision. =#
	## faster arrays
	x,  y, vx, vy, r0, N  =  dat.xy[:,1],  dat.xy[:,2], dat.vxy[:,1], dat.vxy[:,2], 1*dat.r0, 1*dat.N
	## parity check
	dat.parity==1 ? parflip = reverse : parflip = identity
	## mulitiparticle collision events
	hasinteracted = Set{Int}([])
	## loop over particle pairs
	for i=parflip(1:N)
		## loop
		for j=parflip(i+1:N)
			## values
			dx  =  x[j] -  x[i]
			dy  =  y[j] -  y[i]
			dvx = vx[j] - vx[i]
			dvy = vy[j] - vy[i]
			## conditions
			isclose = sqrt(dx^2 + dy^2) < 2*r0
			inbounds = 0<x[i]<1 && 0<y[i]<1 && 0<x[j]<1 && 0<y[j]<1
			## pass if sufficiently close
			if isclose && inbounds
				collide!(dat, i, j)
				i in hasinteracted || j in hasinteracted ? dat.mcc+=1 : nothing
				push!(hasinteracted, i, j)
			end
		end
	end
	## change parity
	dat.parity = (dat.parity + 1) % 2
end

	