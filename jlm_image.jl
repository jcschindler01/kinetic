"""
To generate system image:

>> julia/jlm

using PackageCompiler

create_sysimage(
	["GLMakie"]; 
	sysimage_path="jlm.so",
	precompile_execution_file="jlm_image.jl",
	)

"""

using GLMakie

f = Figure()
a = Axis(f[1,1])

x = 1:10

l = lines!(a, x)
l = lines!(a, x, x.^2)

f

