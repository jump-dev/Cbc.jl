
module Cbc

using BinDeps
@BinDeps.load_dependencies



include("CoinMPInterface.jl")
include("CbcSolverInterface.jl")

using Cbc.CbcMathProgSolverInterface
export CbcSolver

end # module
