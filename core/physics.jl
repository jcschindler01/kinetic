

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
