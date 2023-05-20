



function evolve!(dat::Datapoint; tprime=nothing, kwargs...)
	#= 
	Evolve dat from current time to new time using dat parameters.
	The kwargs can be used to update fields before evolving.
	The tprime can be used to set new time to evolve to, it overrides dt.
	dt = total time to evolve by
	div = number internal substeps
	DT = dt/div = internal timestep
	=#
	## kwargs
	#dump(kwargs)
	## note origin
	dat.origin = "evolve"
	## timestep
	if tprime isa Real
		dat.dt = tprime - dat.t
	end
	## get stepper
	step! = getfield(Kinetic, Symbol("$(dat.integrator)_step!"))
	## internal time substep
	DT = dat.dt/dat.div
	## take steps
	for divs = 1:dat.div
		step!(dat, DT)
	end
	## return
	return dat
end

function free_step!(dat::Datapoint, dt)
	#= Step dat from time [t] to time [t + dt]. =#
	free!(dat, dt)
	wallz!(dat)
end

function naive_step!(dat::Datapoint, dt)
	#= Step dat from time [t] to time [t + dt]. =#
	free!(dat, dt)
	wallz!(dat)
	naivecol!(dat)
end

function free!(dat::Datapoint, dt)
	## free evolution for time dt
	dat.t    += dt
	dat.xy  .+= dt .* dat.vxy
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

function naivecol!(dat::Datapoint)
	#= Collide if close and approaching. (Much faster in current form.) =#
	## faster arrays
	x,  y, vx, vy, r0  =  dat.xy[:,1],  dat.xy[:,2], dat.vxy[:,1], dat.vxy[:,2], 1*dat.r0
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
			## pass if sufficiently close
			if isclose && isapproaching
				collide!(dat, i, j)
			end
		end
	end
end




