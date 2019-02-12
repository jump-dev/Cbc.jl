module Cbc

# This 'using' is required to suppress a warning about Clp not having Libdl in
# its dependencies (Libdl is used by BinaryProvider), e.g.:
# bicycle1885/CodecZlib.jl#26.
using Libdl

if isfile(joinpath(@__DIR__, "..", "deps", "deps.jl"))
    include(joinpath(@__DIR__, "..", "deps", "deps.jl"))
else
    error("Cbc not properly installed. Please run Pkg.build(\"Cbc\")")
end

include("CbcCInterface.jl")
include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

using Cbc.CbcMathProgSolverInterface
export CbcSolver

end
