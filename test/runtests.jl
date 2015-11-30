using Cbc

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
include(joinpath(dirname(@__FILE__),"..","examples","knapsack.jl"))

mixintprogtest(CbcSolver(logLevel=1))

include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(CbcSolver())

if isdir(Pkg.dir("JuMP"))
    include(joinpath(Pkg.dir("JuMP"),"test","runtests.jl"))
end
