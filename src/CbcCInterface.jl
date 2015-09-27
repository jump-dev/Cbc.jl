module CbcCInterface

using Compat

import ..Cbc

export CbcModel,
    loadProblem,
    readMps,
    writeMps,
    setInitialSolution,
    problemName,
    setProblemName,
    getNumElements,
    getVectorStarts,
    getIndices,
    getElements,
    getNumRows,
    getNumCols,
    setObjSense,
    getObjSense,
    getRowLower,
    getRowUpper,
    getRowActivity,
    getColLower,
    getColUpper,
    getObjCoefficients,
    getColSolution,
    setRowUpper,
    setRowLower,
    setObjCoeff,
    setColLower,
    setColUpper,
    isInteger,
    setContinuous,
    setInteger,
    setParameter,
    solve,
    sumPrimalInfeasibilities,
    numberPrimalInfeasibilities,
    getIterationCount,
    isAbandoned,
    isProvenOptimal,
    isProvenInfeasible,
    isContinuousUnbounded,
    isNodeLimitReached,
    isSecondsLimitReached,
    isSolutionLimitReached,
    isInitialSolveAbandoned,
    isInitialSolveProvenOptimal,
    isInitialSolveProvenPrimalInfeasible,
    getObjValue,
    getBestPossibleObjValue,
    getNodecount,
    status,
    secondaryStatus,
    addSOS

    

# helper macros/functions

macro cbc_ccall(func, args...)
    f = "Cbc_$(func)"
    quote
        ccall(($f,Cbc.libcbcsolver), $(args...))
    end
end

type CbcModel
    p::Ptr{Void}
    function CbcModel()
        p = @cbc_ccall newModel Ptr{Void} ()
        prob = new(p)
        finalizer(prob, deleteModel)
        return prob
    end
end

function deleteModel(prob::CbcModel)
    if prob.p == C_NULL
        return
    end
    @cbc_ccall deleteModel Void (Ptr{Void},) prob.p
    prob.p = C_NULL
    return
end

function check_problem(prob::CbcModel)
   if prob.p == C_NULL
       error("Invalid CbcModel")
   end
   return true
end

macro getproperty(T, name)
    @eval function ($name)(prob::CbcModel)
        check_problem(prob)
        @cbc_ccall $name $T (Ptr{Void},) prob.p
    end
end

# Note: we assume COIN_BIG_INDEX and COIN_BIG_DOUBLE
# were not defined when compiling Cbc (which is the
# default)
typealias CoinBigIndex Int32
typealias CoinBigDouble Float64

# copied from Clp.jl
@compat typealias VecOrNothing Union{Vector,Void}
function vec_or_null{T}(::Type{T}, a::VecOrNothing, len::Integer)
    if isequal(a, nothing) || length(a) == 0
        return C_NULL
    else # todo: helpful message if convert fails
        if length(a) != len
            error("Expected vector to have length $len")
        end
        return convert(Vector{T},a)
    end
end

function getVersion()
    s = @cbc_ccall getVersion Ptr{UInt8} ()
    return bytestring(s)
end

function loadProblem(prob::CbcModel,
    constraint_matrix::AbstractMatrix,
    col_lb::VecOrNothing, col_ub::VecOrNothing,
    obj::VecOrNothing,
    row_lb::VecOrNothing, row_ub::VecOrNothing)
    check_problem(prob)

    mat = convert(SparseMatrixCSC{Float64,Int32},constraint_matrix)
    nrow,ncol = size(mat)

    @cbc_ccall loadProblem Void (Ptr{Void}, Int32, Int32, Ptr{CoinBigIndex},
        Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}) prob.p ncol nrow mat.colptr.-@compat(Int32(1)) mat.rowval.-@compat(Int32(1)) mat.nzval vec_or_null(Float64, col_lb, ncol) vec_or_null(Float64, col_ub, ncol) vec_or_null(Float64, obj, ncol) vec_or_null(Float64, row_lb, nrow) vec_or_null(Float64, row_ub, nrow)
end

function readMps(prob::CbcModel, filename::ASCIIString)
    check_problem(prob)
    @cbc_ccall readMps Cint (Ptr{Void}, Ptr{UInt8}) prob.p filename
end

function writeMps(prob::CbcModel, filename::ASCIIString)
    check_problem(prob)
    @cbc_ccall writeMps Cint (Ptr{Void}, Ptr{UInt8}) prob.p filename
end

function setInitialSolution(prob::CbcModel, array::Vector{Float64})
    check_problem(prob)
    @cbc_ccall setInitialSolution Void (Ptr{Void}, Ptr{Float64}) prob.p array
end

function problemName(prob::CbcModel)
    check_problem(prob)
    a = Array(UInt8, 100)
    @cbc_ccall problemName Void (Ptr{Void},Cint,Ptr{UInt8}) prob.p 100, a
    return bytestring(a)
end

function setProblemName(prob::CbcModel)
    check_problem(prob)
    @cbc_ccall setProblemName Cint (Ptr{Void},Ptr{UInt8}) prob.p bytestring(a)
end

function getNumElements(prob::CbcModel)
    check_problem(prob)
    @cbc_ccall getNumElements Cint (Ptr{Void},) prob.p
end

function getVectorStarts(prob::CbcModel)
    check_problem(prob)
    p = @cbc_ccall getVectorStarts Ptr{CoinBigIndex} (Ptr{Void},) prob.p
    num_cols = @compat(Int(getNumCols(prob)))
    return copy(pointer_to_array(p,(num_cols+1,)))
end

function getIndices(prob::CbcModel)
    check_problem(prob)
    p = @cbc_ccall getIndices Ptr{Cint} (Ptr{Void},) prob.p
    nnz = @compat(Int(getNumElements(prob)))
    return copy(pointer_to_array(p,(nnz,)))
end

function getElements(prob::CbcModel)
    check_problem(prob)
    p = @cbc_ccall getElements Ptr{Float64} (Ptr{Void},) prob.p
    nnz = @compat(Int(getNumElements(prob)))
    return copy(pointer_to_array(p,(nnz,)))
end

# TODO: maxNameLength getRowName getColName setColName setRowName

@getproperty Cint getNumRows
@getproperty Cint getNumCols


# 1 : minimize, -1 : maximize
function setObjSense(prob::CbcModel, sense)
    check_problem(prob)
    @cbc_ccall setObjSense Void (Ptr{Void}, Float64) prob.p sense
end

@getproperty Float64 getObjSense

for s in (:getRowLower, :getRowUpper, :getRowActivity)
    @eval function ($s)(prob::CbcModel)
        check_problem(prob)
        nrow = @compat(Int(getNumRows(prob)))
        p = @cbc_ccall $s Ptr{Float64} (Ptr{Void},) prob.p
        return copy(pointer_to_array(p,(nrow,)))
    end
end

for s in (:getColLower, :getColUpper, :getObjCoefficients, :getColSolution)
    @eval function ($s)(prob::CbcModel)
        check_problem(prob)
        ncol = @compat(Int(getNumCols(prob)))
        p = @cbc_ccall $s Ptr{Float64} (Ptr{Void},) prob.p
        return copy(pointer_to_array(p,(ncol,)))
    end
end

for s in (:setRowUpper, :setRowLower, :setObjCoeff, :setColLower, :setColUpper)
    @eval function($s)(prob::CbcModel, index::Integer, value::Float64)
        check_problem(prob)
        @cbc_ccall $s Void (Ptr{Void}, Cint, Float64) prob.p index value
    end
end

function isInteger(prob::CbcModel, index::Integer)
    check_problem(prob)
    v = @cbc_ccall isInteger Cint (Ptr{Void},Cint) prob.p index
    return v == 1
end

function setContinuous(prob::CbcModel, index::Integer)
    check_problem(prob)
    @cbc_ccall setContinuous Void (Ptr{Void},Cint) prob.p index
end

function setInteger(prob::CbcModel, index::Integer)
    check_problem(prob)
    @cbc_ccall setInteger Void (Ptr{Void},Cint) prob.p index
end

function Base.copy(prob::CbcModel)
    p = @cbc_ccall clone Ptr{Void} (Ptr{Void},) prob.p
    prob = CbcModel(p)
    finalizer(prob, deleteModel)
    return prob
end

function setParameter(prob::CbcModel, name::AbstractString, value::AbstractString)
    @cbc_ccall setParameter Void (Ptr{Void},Ptr{UInt8},Ptr{UInt8}) prob.p bytestring(name) bytestring(value)
end

# TODO: registerCallBack clearCallBack

function solve(prob::CbcModel)
    @cbc_ccall solve Cint (Ptr{Void},) prob.p
end

@getproperty Float64 sumPrimalInfeasibilities
@getproperty Cint numberPrimalInfeasibilities

# TODO: checkSolution

@getproperty Cint getIterationCount

for s in (:isAbandoned, :isProvenOptimal, :isProvenInfeasible, :isContinuousUnbounded,
    :isNodeLimitReached, :isSecondsLimitReached, :isSolutionLimitReached, :isInitialSolveAbandoned, :isInitialSolveProvenOptimal, :isInitialSolveProvenPrimalInfeasible)
    @eval function ($s)(prob::CbcModel)
        check_problem(prob)
        v = @cbc_ccall $s Cint (Ptr{Void},) prob.p
        return v != 0
    end
end

@getproperty Float64 getObjValue
@getproperty Float64 getBestPossibleObjValue
@getproperty Cint getNodeCount

# TODO: printSolution

function addSOS(prob::CbcModel, numRows::Integer, rowStarts::Vector{Cint},
    colIndices::Vector{Cint}, weights::Vector{Float64}, typ::Integer)

    @cbc_ccall(addSOS,Void,(Ptr{Void}, Cint, Ptr{Cint}, Ptr{Cint},
                            Ptr{Float64}, Cint),prob.p,numRows,
                            rowStarts-convert(Cint,1),
                            colIndices-convert(Cint,1),
                            weights, typ)

end

# see Cbc_C_Interface.h documentation
@getproperty Cint status
@getproperty Cint secondaryStatus

end
