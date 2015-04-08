using Cbc

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"example","knapsack.jl"))

mixintprogtest(CbcSolver(logLevel=1))

if isdir(Pkg.dir("JuMP"))
    include(joinpath(Pkg.dir("JuMP"),"test","runtests.jl"))
end
