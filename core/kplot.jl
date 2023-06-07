


module KPlot

export Boxplot, update!

using GLMakie
using StatsBase: Histogram, fit, normalize

const vedges = range(0,3,51)
const vbins = (vedges[1:end-1]+vedges[2:end])/2

const boxtheme = Theme(
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
    Axis(fig[2,2],
    )
    return fig
end

@kwdef mutable struct Boxplot
    fig::Figure
    t::Observable{Real} = Observable(0.0)
    x::Observable{Vector{Real}}  = Observable([])
    y::Observable{Vector{Real}}  = Observable([])
    vx::Observable{Vector{Real}} = Observable([])
    vy::Observable{Vector{Real}} = Observable([])
    v::Observable{Vector{Real}}  = Observable([])
    vhist::Observable{Vector{Real}} = 0*vbins
    ann::Observable{String} = "t="
end

function Boxplot()
    bp = Boxplot(fig=newfig())
    scatter!(bp.fig.content[1], bp.x,  bp.y)
    scatter!(bp.fig.content[2], bp.vx, bp.vy)
    scatter!(bp.fig.content[3], vbins, bp.vhist,
        markersize = 10,
        color = :darkslategray3,
        )
    text!(bp.fig.content[4], .5,.5; text=bp.ann, align=(:center,:center))
    return bp
end

function update!(bp::Boxplot, dat)
    bp.t[] = dat.t
    bp.x.val  =  dat.xy[:,1]
    bp.y[]    =  dat.xy[:,2]
    bp.vx.val = dat.vxy[:,1]
    bp.vy[]   = dat.vxy[:,2]
    bp.v.val  = sqrt.(bp.vx.val.^2 + bp.vy.val.^2)
    bp.vhist[] = normalize(fit(Histogram, bp.v.val, vedges); mode=:pdf).weights
    bp.ann[] = "t=$(round(dat.t; digits=3))"
    return bp
end



end




