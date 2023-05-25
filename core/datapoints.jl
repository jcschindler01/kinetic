
@kwdef mutable struct Datapoint
	## parameters
	origin::String="init"					## how it was generated
	N::Int=3 								## number particles
	r0::Float64=0.0							## ball radius (0 = no interaction)
	integrator::String="free"				## collision physics used to generate it
	ic::String="random"						## initial condition it originated from
	dt::Float64=0.1							## evolution time used to generate it
	div::Int=1								## internal timesteps used to generate it
	## counters
	cc::Int=0 								## collision counter
	mcc::Int=0 								## multiparticle-collision counter
	## state								
	t::Float64=0.0							## time
	parity::Int=0							## parity flag
	xy::Array{Float64}=getic(ic,N)[1]		## position
	vxy::Array{Float64}=getic(ic,N)[2]		## velocity
end

