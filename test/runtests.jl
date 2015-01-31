using Cbc

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))

mixintprogtest(CbcSolver(logLevel=1))

if isdir(Pkg.dir("JuMP"))
    include(joinpath(Pkg.dir("JuMP"),"test","runtests.jl"))
end
