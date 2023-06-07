"""
To generate system image:

>> julia

using PackageCompiler

create_sysimage(
	["GLMakie"]; 
	sysimage_path="jlm.so",
	precompile_execution_file="jlm_image.jl",
	)

"""

using GLMakie

boxtheme = Theme(
    fontsize=14,
    palette = (color=[:black,],),
    Axis = (
        limits = (0,1,0,1),
        xticks = [0,1],
        yticks = [0,1],
        aspect = 1,
        xgridvisible = false,
        ygridvisible = false,
    ),
    Lines = (
        color = :black,
    ),
    Scatter = (
        marker = :circle,
        markersize = 5,
    ),
)

mutable struct Boxplot
    fig::Figure
    t::Observable{Real}
    x::Observable{Vector{Real}}
    y::Observable{Vector{Real}}
    vx::Observable{Vector{Real}}
    vy::Observable{Vector{Real}}
end

function Boxplot()
    bp = Boxplot(newfig(), Observable(0.0), Observable([]), Observable([]), Observable([]), Observable([]))
    scatter!(bp.fig.content[1], bp.x,  bp.y)
    scatter!(bp.fig.content[2], bp.vx, bp.vy)
    return bp
end

function newfig()
    ## theme
    set_theme!(boxtheme)
    ## fig
    fig = Figure(resolution = (800,800))
    ## axes
    Axis(fig[1,1], 
        xlabel=L"x", 
        ylabel=L"y",
    )
    Axis(fig[1,2], 
        xlabel=L"v_x", 
        ylabel=L"v_y", 
        limits=(-2,2,-2,2), 
        xticks=-2:2,
        yticks=-2:2,
    )
    Axis(fig[2,1], 
        xlabel=L"v", 
        ylabel=L"p(v)", 
        limits=(0,3,0,2), 
        xticks=0:3,
        yticks=0:2,
    )
    return fig
end


set_theme!(boxtheme)
fig = newfig()
bp = Boxplot()

display(bp.fig)
sleep(2)

