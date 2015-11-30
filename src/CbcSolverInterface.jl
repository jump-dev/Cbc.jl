module CbcMathProgSolverInterface

using Cbc.CbcCInterface

importall MathProgBase.SolverInterface


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
    getvartype,
    setvartype!,
    optimize!,
    status,
    getobjval,
    getobjgap,
    getobjbound,
    getsolution,
    getrawsolver


type CbcMathProgModel <: AbstractLinearQuadraticModel
    inner::CbcModel
    binaries::Vector{Int} # indices of binary variables
end

immutable CbcSolver <: AbstractMathProgSolver
    options 
end
CbcSolver(;kwargs...) = CbcSolver(kwargs)


function CbcMathProgModel(;options...)
    c = CbcModel()
    setParameter(c, "log", "0")
    old_parameters = [:MipMaxSeconds, :LogLevel, :MipMaxSolutions, :MipMaxNodes, :MipAllowableGap, :MipFractionalGap]
    for (optname, optval) in options
        if optname in old_parameters
            warn("Option $optname is no longer recognized. See https://github.com/JuliaOpt/Cbc.jl for renamed list of options.")
        else
            setParameter(c, string(optname), string(optval))
        end
    end
    return CbcMathProgModel(c, Int[])
end

LinearQuadraticModel(s::CbcSolver) = CbcMathProgModel(;s.options...)

ConicModel(s::CbcSolver) = LPQPtoConicBridge(LinearQuadraticModel(s))
supportedcones(s::CbcSolver) = [:Free,:Zero,:NonNeg,:NonPos]

function loadproblem!(m::CbcMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    loadProblem(m.inner, A, collb, colub, obj, rowlb, rowub)
    setsense!(m, sense)
end

function writeproblem(m::CbcMathProgModel, filename::AbstractString)
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

function getvartype(m::CbcMathProgModel)
    ncol = numvar(m)
    vartype = fill(:Cont,ncol)
    for i in 1:ncol
        if isInteger(m.inner, i-1)
            vartype[i] = :Int
        end
    end
    for k in m.binaries
        vartype[k] = :Bin
    end
    return vartype
end

function setvartype!(m::CbcMathProgModel,vartype::Vector{Symbol})
    ncol = numvar(m)
    @assert length(vartype) == ncol
    m.binaries = Int[]
    for i in 1:ncol
        if vartype[i] == :Int
            setInteger(m.inner, i-1)
        elseif vartype[i] == :Cont
            setContinuous(m.inner, i-1)
        elseif vartype[i] == :Bin
            setInteger(m.inner, i-1)
            push!(m.binaries, i)
        else
            error("Unrecognized variable type $(vartype[i])")
        end
    end
end

function optimize!(m::CbcMathProgModel)
    lb = getColLower(m.inner)
    ub = getColUpper(m.inner)
    # tighten bounds on binaries
    for i in m.binaries
        setColLower(m.inner, i-1, max(lb[i],0.0))
        setColUpper(m.inner, i-1, min(ub[i],1.0))
    end
    solve(m.inner)
end

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

function getconstrmatrix(m::CbcMathProgModel)
    starts = getVectorStarts(m.inner)
    rowval = getIndices(m.inner)
    nzval = getElements(m.inner)
    return SparseMatrixCSC(numconstr(m), numvar(m), starts+1, rowval+1, nzval)
end

getobjval(m::CbcMathProgModel) = getObjValue(m.inner)

getobjbound(m::CbcMathProgModel) = getBestPossibleObjValue(m.inner)

function getobjgap(m::CbcMathProgModel)
    b = getBestPossibleObjValue(m.inner)
    f = getObjValue(m.inner)
    return abs(b-f)/abs(f)
end

getsolution(m::CbcMathProgModel) = getColSolution(m.inner)

getrawsolver(m::CbcMathProgModel) = m.inner

function setwarmstart!(m::CbcMathProgModel, v)
    if any(isnan, v)
        Base.warn_once("Ignoring partial starting solution. Cbc requires a feasible value to be specified for all variables.")
        return
    end

    # ignore if not feasible
    @assert length(v) == numvar(m)
    l = getColLower(m.inner)
    u = getColUpper(m.inner)
    for i in 1:length(v)
        if !(l[i] - 1e-6 <= v[i] <= u[i] + 1e-6)
            return
        end
        if isInteger(m.inner, i-1) && !isinteger(l[i])
            return
        end
    end
    lb = getRowLower(m.inner)
    ub = getRowUpper(m.inner)
    A = getconstrmatrix(m)
    rowval = A*v
    for i in 1:numconstr(m)
        if !(lb[i] - 1e-6 <= rowval[i] <= ub[i] + 1e-6)
            return
        end
    end
    setInitialSolution(m.inner, v)
end

function addsos1!(m::CbcMathProgModel, idx, weight::Vector{Float64})
    idxc = convert(Vector{Cint}, idx)
    addSOS(m.inner, 1, Cint[1,length(idx)+1], idxc, weight, 1)
end

function addsos2!(m::CbcMathProgModel, idx, weight::Vector{Float64})
    idxc = convert(Vector{Cint}, idx)
    addSOS(m.inner, 1, Cint[1,length(idx)+1], idxc, weight, 2)
end

end
