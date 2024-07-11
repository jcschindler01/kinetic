
prefix = "/home/jcs/code/kinetic/"

include(prefix * "core/kinetic.jl")
include(prefix * "core/kplot.jl")

using GLMakie
using .Kinetic
using .KPlot


f1 = "runs/txt/hotcold_naive_N500_T5_1720045485.txt"
f2 = "runs/txt/corner_naive_N500_T5_1720110742.txt"
f3 = "runs/txt/gun_naive_N500_T5_1720111984.txt"
f4 = "runs/txt/chain_naive_N500_T5_1720109181.txt"


free1 = "runs/txt/hotcold_free_N500_T5_1720426448.txt"
free2 = "runs/txt/corner_free_N500_T5_1720111905.txt"
free3 = "runs/txt/gun_free_N500_T5_1720112713.txt"
free4 = "runs/txt/chain_free_N500_T5_1720109831.txt"

runs = [f1,f2,f3,f4]
free = [free1,free2,free3,free4]

cd(@__DIR__)

for data in [f2]
	record_animation(prefix * data, 
		save=false, 
		fps=20, maxlines=100, delay=25,
		fmt="gif"
		)
end


