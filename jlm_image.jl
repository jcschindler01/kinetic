"""
To generate system image:

>> julia/jlm

using PackageCompiler

create_sysimage(
	["GLMakie"]; 
	sysimage_path="jlm.so",
	precompile_execution_file="jlm_image.jl",
	)

## First time will use tons of RAM (~16+GB) and will crash if RAM runs out.
## Add swap!!
## Will also take lots of time (~1.5+hrs).
## Set keep awake!!

"""


include("runs/run.jl")

