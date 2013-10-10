
export CoinMPSolver,
    model,
    loadproblem,
    writeproblem,
    updatemodel,
    setsense,
    getsense,
    numvar,
    numconstr,
    setvartype,
    optimize,
    status,
    getobjval,
    getobjbound,
    getsolution,
    getrawsolver


type CoinMPSolver <: LinprogSolver
    inner::CoinProblem
end

function model(;options...)
    c = CoinProblem()
    setOption(c, "LogLevel", 0)
    for (optname, optval) in options
        setOption(c, string(optname), optval)
    end
    return CoinMPSolver(c)
end

loadproblem(m::CoinMPSolver, A, collb, colub, obj, rowlb, rowub) =
    LoadMatrix(m.inner, 1, 0.0, obj, collb, colub, rowlb, rowub, A)

function writeproblem(m::CoinMPSolver, filename::String)
    if endswith(filename,".mps")
        WriteFile(m.inner, 3, filename)
    else
        error("Only MPS output supported")
    end
end

updatemodel(m::CoinMPSolver) = nothing

function setsense(m::CoinMPSolver,sense)
    if sense != :Min
        error("Only minimization sense currently supported")
    end
end

getsense(m::CoinMPSolver) = :Min

numvar(m::CoinMPSolver) = GetColCount(m.inner)
numconstr(m::CoinMPSolver) = GetRowCount(m.inner)

function setvartype(m::CoinMPSolver,vartype)
    ncol = numvar(m)
    @assert length(vartype) == ncol
    coltype = Array(Uint8,ncol)
    for i in 1:ncol
        @assert vartype[i] == 'I' || vartype[i] == 'C'
        coltype[i] = vartype[i]
    end
    LoadInteger(m.inner,coltype)
end

optimize(m::CoinMPSolver) = OptimizeProblem(m.inner)

function status(m::CoinMPSolver)
    stat = GetSolutionText(m.inner)
    # CoinMP status reporting is faulty,
    # add logic from CbcModel.cpp.
    objval = GetObjectValue(m.inner)
    if stat == "Optimal solution found" 
        if objval < 1e30
            return :Optimal
        else
            return :Infeasible
        end
    elseif stat == "Problem primal infeasible"
        return :Infeasible
    elseif stat == "Problem dual infeasible" # what does this mean for MIP??
        return :Unbounded
    elseif stat == "Stopped on iterations" || stat == "Stopped by user"
        return :UserLimit
    elseif stat == "Stopped due to errors"
        return :Error
    else
        error("Internal library error")
    end
end

getobjval(m::CoinMPSolver) = GetObjectValue(m.inner)

getobjbound(m::CoinMPSolver) = GetMipBestBound(m.inner)

getsolution(m::CoinMPSolver) = GetSolutionValues(m.inner)

getrawsolver(m::CoinMPSolver) = m.inner
