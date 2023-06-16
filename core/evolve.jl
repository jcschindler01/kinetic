



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
	free!(dat, dt/2)
	wallv!(dat)
	free!(dat, dt/2)
end

function naive_step!(dat::Datapoint, dt)
	#= Step dat from time [t] to time [t + dt]. =#
	free!(dat, dt/2)
	wallv!(dat)
	naivecol!(dat)
	free!(dat, dt/2)
end

function sym_step!(dat::Datapoint, dt)
	#= Step dat from time [t] to time [t + dt]. =#
	free!(dat, dt/2)
	wallv!(dat)
	symcol!(dat)
	free!(dat, dt/2)
end

function sym1_step!(dat::Datapoint, dt)
	#= Step dat from time [t] to time [t + dt]. =#
	free!(dat, dt/2)
	wallv!(dat)
	symcol1!(dat)
	free!(dat, dt/2)
end
