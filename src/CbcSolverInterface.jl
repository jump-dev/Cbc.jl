module CbcMathProgSolverInterface

using Cbc.CbcCInterface

require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
importall MathProgSolverInterface


export CbcMathProgModel,
    CbcSolver,
    model,
    loadproblem!,
    writeproblem,
    updatemodel!,
    setsense!,
    getsense,
    numvar,
    numconstr,
    setvartype!,
    optimize!,
    status,
    getobjval,
    getobjbound,
    getsolution,
    getrawsolver


type CbcMathProgModel <: AbstractMathProgModel
    inner::CbcModel
end

immutable CbcSolver <: AbstractMathProgSolver
    options 
end
CbcSolver(;kwargs...) = CbcSolver(kwargs)


function CbcMathProgModel(;options...)
    c = CbcModel()
    setParameter(c, "log", "0")
    for (optname, optval) in options
        setParameter(c, string(optname), string(optval))
    end
    return CbcMathProgModel(c)
end

model(s::CbcSolver) = CbcMathProgModel(;s.options...)

function loadproblem!(m::CbcMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    loadProblem(m.inner, A, collb, colub, obj, rowlb, rowub)
    setsense!(m, sense)
end

function writeproblem(m::CbcMathProgModel, filename::String)
    if endswith(filename,".mps")
        writeMps(m.inner, filename)
    else
        error("Only MPS output supported")
    end
end

updatemodel(m::CbcMathProgModel) = nothing

function setsense!(m::CbcMathProgModel,sense)
    @assert sense == :Min || sense == :Max
    if sense == :Min
       setObjSense(m.inner, 1)
    else
       setObjSense(m.inner, -1)
    end
end

function getsense(m::CbcMathProgModel)
    s = getObjSense(m.inner)
    if s == 1
        return :Min
    elseif s == -1
        return :Max
    else
        error("Internal error: Unknown sense $s")
    end
end

numvar(m::CbcMathProgModel) = getNumCols(m.inner)
numconstr(m::CbcMathProgModel) = getNumRows(m.inner)

function setvartype!(m::CbcMathProgModel,vartype)
    ncol = numvar(m)
    @assert length(vartype) == ncol
    for i in 1:ncol
        @assert vartype[i] == 'I' || vartype[i] == 'C'
        if vartype[i] == 'I'
            setInteger(m.inner, i-1)
        else
            setContinuous(m.inner, i-1)
        end
    end
end

optimize!(m::CbcMathProgModel) = solve(m.inner)

function status(m::CbcMathProgModel)
    if isProvenOptimal(m.inner)
        return :Optimal
    elseif isProvenInfeasible(m.inner)
        return :Infeasible
    elseif isContinuousUnbounded(m.inner)
        return :Unbounded # is this correct?
    elseif isNodeLimitReached(m.inner) || isSecondsLimitReached(m.inner) || isSolutionLimitReached(m.inner)
        return :UserLimit
    elseif isAbandoned(m.inner)
        return :Error
    else
        error("Internal error: Unrecognized solution status")
    end
end

getobjval(m::CbcMathProgModel) = getObjValue(m.inner)

getobjbound(m::CbcMathProgModel) = getBestPossibleObjValue(m.inner)

getsolution(m::CbcMathProgModel) = getColSolution(m.inner)

getrawsolver(m::CbcMathProgModel) = m.inner

setwarmstart!(m::CbcMathProgModel, v) = setInitialSolution(m.inner, v)

end
