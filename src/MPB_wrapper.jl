module CbcMathProgSolverInterface

using Cbc.CbcCInterface

using SparseArrays

import MathProgBase
const MPB = MathProgBase

export CbcMathProgModel, CbcSolver

mutable struct CbcMathProgModel <: MPB.AbstractLinearQuadraticModel
    inner::CbcModel
    check_warmstart::Bool
    binaries::Vector{Int} # indices of binary variables
end

mutable struct CbcSolver <: MPB.AbstractMathProgSolver
    options
end
CbcSolver(;kwargs...) = CbcSolver(kwargs)


function CbcMathProgModel(;check_warmstart::Bool=true,options...)
    c = CbcModel()
    setParameter(c, "logLevel", "0")
    for (optname, optval) in options
        setParameter(c, string(optname), string(optval))
    end
    return CbcMathProgModel(c, check_warmstart, Int[])
end

MPB.LinearQuadraticModel(s::CbcSolver) = CbcMathProgModel(;s.options...)

MPB.ConicModel(s::CbcSolver) = MPB.LPQPtoConicBridge(MPB.LinearQuadraticModel(s))
MPB.supportedcones(s::CbcSolver) = [:Free,:Zero,:NonNeg,:NonPos]

function MPB.setparameters!(s::CbcSolver; mpboptions...)
    opts = collect(Any, s.options)
    silent = false
    for (optname, optval) in mpboptions
        if optname == :TimeLimit
            @static if VERSION >= v"0.7-"
                push!(opts, :seconds => optval)
            else
                push!(opts, (:seconds, optval))
            end
        elseif optname == :Silent
            if optval == true
                @static if VERSION >= v"0.7-"
                    push!(opts, :logLevel => 0)
                else
                    push!(opts, (:logLevel, 0))
                end
            end
        else
            error("Unrecognized parameter $optname")
        end
    end
    s.options = opts
end

function MPB.setparameters!(s::CbcMathProgModel; mpboptions...)
    for (optname, optval) in mpboptions
        if optname == :TimeLimit
            setParameter(s.inner, "seconds", string(optval))
        elseif optname == :Silent
            if optval == true
                setParameter(s.inner, "logLevel","0")
            end
        else
            error("Unrecognized parameter $optname")
        end
    end
end

function MPB.loadproblem!(m::CbcMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    loadProblem(m.inner, A, collb, colub, obj, rowlb, rowub)
    MPB.setsense!(m, sense)
end

function MPB.writeproblem(m::CbcMathProgModel, filename::AbstractString)
    if endswith(filename,".mps")
        writeMps(m.inner, filename)
    else
        error("Only MPS output supported")
    end
end

function MPB.setsense!(m::CbcMathProgModel,sense)
    @assert sense == :Min || sense == :Max
    if sense == :Min
       setObjSense(m.inner, 1)
    else
       setObjSense(m.inner, -1)
    end
end

function MPB.getsense(m::CbcMathProgModel)
    s = getObjSense(m.inner)
    if s == 1
        return :Min
    elseif s == -1
        return :Max
    else
        error("Internal error: Unknown sense $s")
    end
end

MPB.numvar(m::CbcMathProgModel) = getNumCols(m.inner)
MPB.numconstr(m::CbcMathProgModel) = getNumRows(m.inner)

function MPB.getvartype(m::CbcMathProgModel)
    ncol = MPB.numvar(m)
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

function MPB.setvartype!(m::CbcMathProgModel,vartype::Vector{Symbol})
    ncol = MPB.numvar(m)
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

function MPB.optimize!(m::CbcMathProgModel)
    lb = getColLower(m.inner)
    ub = getColUpper(m.inner)
    # tighten bounds on binaries
    for i in m.binaries
        setColLower(m.inner, i-1, max(lb[i],0.0))
        setColUpper(m.inner, i-1, min(ub[i],1.0))
    end
    solve(m.inner)
end

function MPB.status(m::CbcMathProgModel)
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

function MPB.getconstrmatrix(m::CbcMathProgModel)
    starts = getVectorStarts(m.inner)
    rowval = getIndices(m.inner)
    nzval = getElements(m.inner)
    return SparseMatrixCSC(MPB.numconstr(m), MPB.numvar(m), starts .+ 1,
                           rowval .+ 1, nzval)
end

MPB.getobjval(m::CbcMathProgModel) = getObjValue(m.inner)

MPB.getobjbound(m::CbcMathProgModel) = getBestPossibleObjValue(m.inner)

function MPB.getobjgap(m::CbcMathProgModel)
    b = getBestPossibleObjValue(m.inner)
    f = getObjValue(m.inner)
    return abs(b-f)/abs(f)
end

MPB.getnodecount(m::CbcMathProgModel) = getNodeCount(m.inner)

MPB.getsolution(m::CbcMathProgModel) = getColSolution(m.inner)

MPB.getrawsolver(m::CbcMathProgModel) = m.inner

function MPB.setwarmstart!(m::CbcMathProgModel, v)
    if any(isnan, v)
        @warn("Ignoring partial starting solution. Cbc requires a feasible " *
              "value to be specified for all variables.", maxlog = 1)
        return
    end

    # ignore if not feasible
    if m.check_warmstart
        @assert length(v) == MPB.numvar(m)
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
        A = MPB.getconstrmatrix(m)
        rowval = A*v
        for i in 1:MPB.numconstr(m)
            if !(lb[i] - 1e-6 <= rowval[i] <= ub[i] + 1e-6)
                return
            end
        end
    end
    setInitialSolution(m.inner, v)
end

function MPB.addsos1!(m::CbcMathProgModel, idx, weight::Vector{Float64})
    idxc = convert(Vector{Cint}, idx)
    addSOS(m.inner, 1, Cint[1,length(idx)+1], idxc, weight, 1)
end

function MPB.addsos2!(m::CbcMathProgModel, idx, weight::Vector{Float64})
    idxc = convert(Vector{Cint}, idx)
    addSOS(m.inner, 1, Cint[1,length(idx)+1], idxc, weight, 2)
end

end
