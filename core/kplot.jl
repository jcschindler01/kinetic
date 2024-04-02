


module KPlot

export Boxplot, update!, showinit!, animate, record_animation

include("kinetic.jl")
using .Kinetic

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
        ylabel=L"$S/N$ (bits)",
        limits=(0,1,0,20),
        xticks=0:40,
        yticks=0:1:20,
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
    S_localE::Observable{Vector{Real}} = Observable([])
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
    hlines!(bp.fig.content[4], bp.S0, label=L"S(\tau)", color = (:black, 0.9), linestyle=:dash)
    lines!(bp.fig.content[4], bp.tt, bp.S_spatial, label=L"M_{\textrm{space}}", color=:darkgreen, alpha=.8)
    lines!(bp.fig.content[4], bp.tt, bp.S_velocity, label=L"M_{\textrm{speed}}", color=:darkblue, alpha=.8)
    lines!(bp.fig.content[4], bp.tt, bp.S_localE, label=L"M_{E_A} \otimes M_{E_B}", color=:cyan, alpha=.8)
    axislegend(bp.fig.content[4], position=:rb, framevisible=false)
    text!(bp.fig.content[3], 2.9, 1.95; text=bp.ann, align=(:right,:top), font="Courier", fontsize=16)
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
    append!(bp.S_localE.val, S_localE(dat)/dat.N)
    bp.tt[] = bp.tt.val
    bp.S_spatial[] = bp.S_spatial.val
    bp.S_velocity[] = bp.S_velocity.val
    bp.S_localE[] = bp.S_localE.val
    limits!(bp.fig.content[4], (nothing, nothing), (nothing, ceil(Int, bp.S0.val)+1/2))
    #autolimits!(bp.fig.content[4])
    #####
    bp.ann[] =  """
                  t=$(round(dat.t; digits=3))

                  N=$(dat.N)
                 r0=$(dat.r0)
                 af=$(round(area(dat); digits=3))
                int=$(dat.integrator)
                 ic=$(dat.ic)

                  E/N=$(round(temp(dat); digits=3))
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


function animate(filename; tmin=-Inf, tmax=Inf, rate=.2, loops=1, delay=1)
    for loop=1:loops
        ## initialize
        dat = Datapoint()
        bp = Boxplot()
        display(bp.fig)
        ## go
        open(filename, "r") do io
            ## header
            readline(io)
            ## initial data
            dat = Datapoint()
            from_qp!(dat, readline(io))
            dat.xy0  = 1 .* dat.xy
            dat.vxy0 = 1 .* dat.vxy
            update!(bp, dat)
            if loop==1
                dump(dat)
                showinit!(bp, dat)
            end
            ## go
            while !eof(io)
                from_qp!(dat, readline(io))
                if (tmin <= dat.t <= tmax)
                    update!(bp,  dat)
                    sleep(dat.dt/rate)
                end
            end
        end
    end
    sleep(delay)
end

function record_animation(filename; fmt="gif", save=false, fps=15, maxlines=300, delay=25)
    ##
    dat = Datapoint()
    bp = Boxplot()
    showinit!(bp, dat)
    ##
    linecount = countlines(filename)-2
    println(linecount)
    linecount = min(linecount, maxlines)
    open(filename, "r") do io
        #### PREP TO RECORD ####
        ## lose header
        readline(io)
        ## initial data
        dat = Datapoint()
        from_qp!(dat, readline(io))
        dat.xy0  = 1 .* dat.xy
        dat.vxy0 = 1 .* dat.vxy
        update!(bp, dat)
        #### READY TO RECORD ###
        ## 
        global tempfile = "temp.$(fmt)"
        global outfile = "$(fmt)/$(dat.ic)_$(dat.integrator)_N$(dat.N)_$(round(Int,time())).$(fmt)"
        ##
        lines = 1:(linecount+delay)
        function step(line)
            println(line)
            if (line>delay) && (!eof(io))
                from_qp!(dat, readline(io))
                update!(bp,  dat)
            end
        end
        record(step, bp.fig, tempfile, lines; framerate=fps)
    end
    if save==true
        cp(tempfile, outfile)
    end
end







end




