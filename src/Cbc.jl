
module Cbc

include("../deps/deps.jl")



include("CoinMPInterface.jl")
include("CbcSolverInterface.jl")

using Cbc.CbcMathProgSolverInterface
export CbcSolver

end # module
