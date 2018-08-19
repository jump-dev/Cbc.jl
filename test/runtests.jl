using Cbc

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
include(joinpath(dirname(@__FILE__),"..","examples","knapsack.jl"))

mixintprogtest(CbcSolver(logLevel=1))

include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
solver = CbcSolver(logLevel=1)
MathProgBase.setparameters!(solver, Silent=true, TimeLimit=100.0)
println("Cbc should not display any output from these tests")
coniclineartest(solver)

include("MOIWrapper.jl")
