
cd(@__DIR__)

include("../core/kinetic.jl")
include("../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot



animate("txt/gun_naive_N500_T4_1707926341.txt", rate=.25)

# animate("txt/chain_naive_N1000_T10_1707918727.txt", rate=.25)

