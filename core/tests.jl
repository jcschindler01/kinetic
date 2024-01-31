
"""
To test.

>>> julia
> include("tests.jl")

You should then be able to edit file and run functions to test.
"""



using Revise

include("kinetic.jl")
include("kplot.jl")
using .Kinetic, .KPlot


## tests module ##
##################
module KTests

export alltests
export test1

include("helpers.jl")


## run all tests
function alltests()
    test1()
end


## test 1: just load everything
function test1()
    print("test1")
end




end
##################
##################


## import module and run tests
using .KTests
alltests()

