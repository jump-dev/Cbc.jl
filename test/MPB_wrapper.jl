@static if VERSION < v"0.7.0-beta2.203"
    include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
else
    import MathProgBase
    include(joinpath(dirname(pathof(MathProgBase)), "..", "test", "mixintprog.jl"))
end
include(joinpath(dirname(@__FILE__),"..","examples","knapsack.jl"))

mixintprogtest(CbcSolver(logLevel=1))

@static if VERSION < v"0.7.0-beta2.203"
    include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
else
    include(joinpath(dirname(pathof(MathProgBase)), "..", "test", "conicinterface.jl"))
end
solver = CbcSolver(logLevel=1)
MathProgBase.setparameters!(solver, Silent=true, TimeLimit=100.0)
println("Cbc should not display any output from these tests")
coniclineartest(solver)
