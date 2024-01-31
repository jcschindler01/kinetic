
# include("../../core/kinetic.jl")
# include("../../core/kplot.jl")

# using GLMakie
# using .Kinetic
# using .KPlot

include("../../core/entropic.jl")

using .KEntropic

N = 100000
xybins = 3
fbins = 50

# x,y =rand(N), rand(N)
# df = 1/fbins

# Nxy = hist2d(x, y, nbins=xybins)
# fxy = Nxy/N
# fxy_int = floor.(Int,fxy./df)
# fxy_coarse = fxy_int.*df
# fxy_string = string(fxy_int)



## monte carlo ##
N = 5000
xybins = 3
fbins = 10
df = 1/fbins

qdict = Dict{String,Int}([("_N",0)])
trials = 100000
for mc in 1:trials
	##
	x,y =rand(N), rand(N)
	Nxy = hist2d(x, y, nbins=xybins)
	fxy = Nxy/N
	fxy_int = floor.(Int,fxy./df)
	fxy_coarse = fxy_int.*df
	fxy_string = string(fxy_int)
	##
	qdict["_N"] += 1
	if fxy_string in keys(qdict)
		qdict[fxy_string] += 1
	else
		qdict[fxy_string]  = 1
	end
end

qdict_sorted = sort(qdict, byvalue=true, rev=true)


for s in collect(keys(qdict_sorted))[2:-1]
	println(s) 
	println(qref(qdict_sorted,s))
	println()
end




