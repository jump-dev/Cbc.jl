using Cbc

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))

mixintprogtest(CbcSolver(LogLevel=1))
