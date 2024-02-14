
cd(@__DIR__)

include("../core/kinetic.jl")
include("../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot



animate("gun_naive_N500_T4_1707926341.txt", tmax=3)


