__precompile__()

module Cbc

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("Cbc not properly installed. Please run Pkg.build(\"Cbc\")")
end



include("CbcCInterface.jl")
include("CbcSolverInterface.jl")

using Cbc.CbcMathProgSolverInterface
export CbcSolver

end # module
