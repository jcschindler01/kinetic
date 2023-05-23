



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
	for kwa in kwargs
		if kwa[1] in fieldnames(Datapoint)
			val = oftype(getfield(dat, kwa[1]), kwa[2])
			setfield!(dat, kwa[1], val)
		end
	end
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

function naivecol!(dat::Datapoint)
	#= Collide if close and approaching. =#
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


function symcol!(dat::Datapoint; parity=0)
	#= Collide if close and approaching. =#
	## faster arrays
	x,  y, vx, vy, r0  =  dat.xy[:,1],  dat.xy[:,2], dat.vxy[:,1], dat.vxy[:,2], 1*dat.r0
	## loop over particle pairs
	for i=1:dat.N
		## parity check
		if parity==1; i=N-i; end
		## loop
		for j=i+1:dat.N
			## parity check
			if parity==1; j=N-j; end
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



