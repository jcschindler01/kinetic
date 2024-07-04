
cd(@__DIR__)

include("../core/kinetic.jl")
include("../core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot



record_animation("txt/hotcold_naive_N500_T5_1720043792.txt", save=true)

# "chain_naive_N1000_T10_1707918727"
# "corner_naive_N500_T4_1707916201"
# "corner_free_N500_T4_1708531699"
# "hotcold_naive_N500_T1_1711556349"
# "gun_naive_N500_T4_1707926341"
# "hotcold_naive_N500_T5_1720043792"