


module KPlot

export Boxplot, update!, showinit!

using GLMakie
using StatsBase: Histogram, fit, normalize

include("helpers.jl")
include("entropic.jl")
maxboltz(v, T::Observable{Real}) = maxboltz(v; T=T.val)

const vedges = range(0,3,31)
const vbins = (vedges[1:end-1]+vedges[2:end])/2
const vsmooth = range(0,3,1001)
const tsmooth = range(0,100,2)

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
        lw = 2,
    ),
    Scatter = (
        marker = :circle,
        markersize = 5,
        strokevisible = false,
    ),
)

function newfig()
    ## theme
    set_theme!(boxtheme)
    ## fig
    fig = Figure(size = (800,800))
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
        xlabel=L"t", 
        ylabel=L"S/N (bits)",
        limits=(0,1,0,20),
        xticks=0:40,
        yticks=0:5:20,
    )
    return fig
end

@kwdef mutable struct Boxplot
    fig::Figure
    t::Observable{Real} = Observable(0.0)
    temp::Observable{Real} = Observable(0.5)
    ms::Observable{Real} = Observable(5)
    x::Observable{Vector{Real}}  = Observable([])
    y::Observable{Vector{Real}}  = Observable([])
    vx::Observable{Vector{Real}} = Observable([])
    vy::Observable{Vector{Real}} = Observable([])
    v::Observable{Vector{Real}}  = Observable([])
    vhist::Observable{Vector{Real}} = 0*vbins
    tt::Observable{Vector{Real}} = Observable([])
    S0::Observable{Real} = Observable(0.0)
    S_spatial::Observable{Vector{Real}} = Observable([])
    S_velocity::Observable{Vector{Real}} = Observable([])
    ann::Observable{String} = "t="
end

function Boxplot()
    bp = Boxplot(fig=newfig())
    scatter!(bp.fig.content[1], bp.x,  bp.y;  markersize=bp.ms)
    scatter!(bp.fig.content[2], bp.vx, bp.vy; markersize=bp.ms)
    scatter!(bp.fig.content[3], vbins, bp.vhist,
        markersize = 10,
        color = :darkslategray3,
        )
    lines!(bp.fig.content[3], vsmooth, maxboltz(vsmooth, bp.temp), -ones(length(vsmooth)),
        color = (:darkslategray3, 0.5),
        )
    hlines!(bp.fig.content[4], bp.S0, color = (:gray, 0.3))
    lines!(bp.fig.content[4], bp.tt, bp.S_spatial, color=:darkgreen)
    lines!(bp.fig.content[4], bp.tt, bp.S_velocity, color=:darkblue)
    #text!(bp.fig.content[4], .2,.6; text=bp.ann, align=(:left,:center), font="Courier", fontsize=16)
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
    bp.temp[] = round(temp(dat); digits=3)
    bp.ms[] = ceil(Int, dotsize(dat.r0))
    ## entropy ##
    bp.S0[] = Stau(N=dat.N)/dat.N
    append!(bp.tt.val, dat.t)
    append!(bp.S_spatial.val, S_spatial(dat)/dat.N)
    append!(bp.S_velocity.val, S_velocity(dat)/dat.N)
    bp.tt[] = bp.tt.val
    bp.S_spatial[] = bp.S_spatial.val
    bp.S_velocity[] = bp.S_velocity.val
    limits!(bp.fig.content[4], (nothing, nothing), (0,20))
    #####
    bp.ann[] =  """
                  t=$(round(dat.t; digits=3))

                  N=$(dat.N)
                 r0=$(dat.r0)
                 af=$(round(area(dat); digits=3))
                int=$(dat.integrator)
                 ic=$(dat.ic)

                  T=$(round(temp(dat); digits=3))
                 cc=$(dat.cc)
                mcc=$(dat.mcc)

                pdn=$(round(phasedist(dat)/dat.N; digits=5))
                """
    return bp
end


function showinit!(bp::Boxplot, dat)
    col = "#dddddd"
    scatter!(bp.fig.content[1],  dat.xy0[:,1],  dat.xy0[:,2], -ones(dat.N); color=col, markersize=bp.ms)
    scatter!(bp.fig.content[2], dat.vxy0[:,1], dat.vxy0[:,2], -ones(dat.N); color=col, markersize=bp.ms)
    v0 = sqrt.(dat.vxy[:,1].^2 + dat.vxy[:,2].^2)
    vhist0 = normalize(fit(Histogram, v0, vedges); mode=:pdf).weights
    scatter!(bp.fig.content[3], vbins, vhist0, -ones(length(vbins)),
        markersize = 10,
        color = col,
        )
end



end




