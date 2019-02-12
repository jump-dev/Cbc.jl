module CbcCInterface

import ..Cbc.libcbcsolver

using SparseArrays

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
    unsafe_getRowLower,
    unsafe_getRowUpper,
    unsafe_getRowActivity,
    getRowLower,
    getRowUpper,
    getRowActivity,
    unsafe_getColLower,
    unsafe_getColUpper,
    unsafe_getObjCoefficients,
    unsafe_getColSolution,
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
    optimizeNotCalled,
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
    getNodeCount,
    status,
    secondaryStatus,
    addSOS



# helper macros/functions

macro cbc_ccall(func, args...)
    args = map(esc, args)
    f = "Cbc_$(func)"
    quote
        ccall(($f, libcbcsolver), $(args...))
    end
end

mutable struct CbcModel
    p::Ptr{Cvoid}
    function CbcModel()
        p = @cbc_ccall newModel Ptr{Cvoid} ()
        prob = new(p)
        finalizer(deleteModel, prob)
        return prob
    end
end

# This is defined so that we can pass a CbcModel as an argument to ccall without worrying about it being garbage collected during the call.
Base.unsafe_convert(::Type{Ptr{Cvoid}}, prob::CbcModel) = prob.p

function deleteModel(prob::CbcModel)
    if prob.p == C_NULL
        return
    end
    @cbc_ccall deleteModel Cvoid (Ptr{Cvoid},) prob
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
        @cbc_ccall $name $T (Ptr{Cvoid},) prob
    end
end

# Note: we assume COIN_BIG_INDEX and COIN_BIG_DOUBLE
# were not defined when compiling Cbc (which is the
# default)
const CoinBigIndex = Int32
const CoinBigDouble = Float64

# copied from Clp.jl
const VecOrNothing = Union{Vector,Cvoid}
function vec_or_null(::Type{T}, a::VecOrNothing, len::Integer) where T
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
    return string(s)
end

function loadProblem(prob::CbcModel,
    constraint_matrix::AbstractMatrix,
    col_lb::VecOrNothing, col_ub::VecOrNothing,
    obj::VecOrNothing,
    row_lb::VecOrNothing, row_ub::VecOrNothing)
    check_problem(prob)

    mat = convert(SparseMatrixCSC{Float64,Int32},constraint_matrix)
    nrow,ncol = size(mat)

    @cbc_ccall loadProblem Cvoid (Ptr{Cvoid}, Int32, Int32, Ptr{CoinBigIndex},
        Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}) prob ncol nrow mat.colptr.-Int32(1) mat.rowval.-Int32(1) mat.nzval vec_or_null(Float64, col_lb, ncol) vec_or_null(Float64, col_ub, ncol) vec_or_null(Float64, obj, ncol) vec_or_null(Float64, row_lb, nrow) vec_or_null(Float64, row_ub, nrow)
end

function readMps(prob::CbcModel, filename::String)
    check_problem(prob)
    @assert isascii(filename)
    @cbc_ccall readMps Cint (Ptr{Cvoid}, Ptr{UInt8}) prob filename
end

function writeMps(prob::CbcModel, filename::String)
    check_problem(prob)
    @assert isascii(filename)
    @cbc_ccall writeMps Cint (Ptr{Cvoid}, Ptr{UInt8}) prob filename
end

function setInitialSolution(prob::CbcModel, array::Vector{Float64})
    check_problem(prob)
    @cbc_ccall setInitialSolution Cvoid (Ptr{Cvoid}, Ptr{Float64}) prob array
end

function problemName(prob::CbcModel)
    check_problem(prob)
    a = Array(UInt8, 100)
    @cbc_ccall problemName Cvoid (Ptr{Cvoid},Cint,Ptr{UInt8}) prob 100 a
    return string(a)
end

#function setProblemName(prob::CbcModel)
#    check_problem(prob)
#    @cbc_ccall setProblemName Cint (Ptr{Cvoid},Ptr{UInt8}) prob bytestring(a)
#end

function getNumElements(prob::CbcModel)
    check_problem(prob)
    @cbc_ccall getNumElements Cint (Ptr{Cvoid},) prob
end

function getVectorStarts(prob::CbcModel)
    check_problem(prob)
    p = @cbc_ccall getVectorStarts Ptr{CoinBigIndex} (Ptr{Cvoid},) prob
    num_cols = Int(getNumCols(prob))
    return copy(unsafe_wrap(Array,p,(num_cols+1,)))
end

function getIndices(prob::CbcModel)
    check_problem(prob)
    p = @cbc_ccall getIndices Ptr{Cint} (Ptr{Cvoid},) prob
    nnz = Int(getNumElements(prob))
    return copy(unsafe_wrap(Array,p,(nnz,)))
end

function getElements(prob::CbcModel)
    check_problem(prob)
    p = @cbc_ccall getElements Ptr{Float64} (Ptr{Cvoid},) prob
    nnz = Int(getNumElements(prob))
    return copy(unsafe_wrap(Array,p,(nnz,)))
end

# TODO: maxNameLength getRowName getColName setColName setRowName

@getproperty Cint getNumRows
@getproperty Cint getNumCols


# 1 : minimize, -1 : maximize
function setObjSense(prob::CbcModel, sense)
    check_problem(prob)
    @cbc_ccall setObjSense Cvoid (Ptr{Cvoid}, Float64) prob sense
end

@getproperty Float64 getObjSense

for s in (:getRowLower, :getRowUpper, :getRowActivity)
    # This macro creates a unsafe wrapper to the following function
    # It is unsafe due to the returned array residing in memory owned by Cbc
    # The unsafe methods should be recalled when ever a change has been made to Cbc
    unsafe_s = Symbol("unsafe", '_', s)
    @eval function ($unsafe_s)(prob::CbcModel)
        check_problem(prob)
        nrow = Int(getNumRows(prob))
        p = @cbc_ccall $s Ptr{Float64} (Ptr{Cvoid},) prob
        return unsafe_wrap(Array, p, (nrow,))
    end
    # By making a copy of the unsafe array we can depend on it not changing.
    @eval ($s)(prob::CbcModel) = copy(($unsafe_s)(prob))
end

for s in (:getColLower, :getColUpper, :getObjCoefficients, :getColSolution)
    # This macro creates a unsafe wrapper to the following function
    # It is unsafe due to the returned array residing in memory owned by Cbc
    # The unsafe methods should be recalled when ever a change has been made to Cbc
    unsafe_s = Symbol("unsafe", '_', s)
    @eval function ($unsafe_s)(prob::CbcModel)
        check_problem(prob)
        ncol = Int(getNumCols(prob))
        p = @cbc_ccall $s Ptr{Float64} (Ptr{Cvoid},) prob
        return unsafe_wrap(Array, p, (ncol,))
    end
    # By making a copy of the unsafe array we can depend on it not changing.
    @eval ($s)(prob::CbcModel) = copy(($unsafe_s)(prob))
end

for s in (:setRowUpper, :setRowLower, :setObjCoeff, :setColLower, :setColUpper)
    @eval function($s)(prob::CbcModel, index::Integer, value::Float64)
        check_problem(prob)
        @cbc_ccall $s Cvoid (Ptr{Cvoid}, Cint, Float64) prob index value
    end
end

function isInteger(prob::CbcModel, index::Integer)
    check_problem(prob)
    v = @cbc_ccall isInteger Cint (Ptr{Cvoid},Cint) prob index
    return v == 1
end

function setContinuous(prob::CbcModel, index::Integer)
    check_problem(prob)
    @cbc_ccall setContinuous Cvoid (Ptr{Cvoid},Cint) prob index
end

function setInteger(prob::CbcModel, index::Integer)
    check_problem(prob)
    @cbc_ccall setInteger Cvoid (Ptr{Cvoid},Cint) prob index
end

function Base.copy(prob::CbcModel)
    p = @cbc_ccall clone Ptr{Cvoid} (Ptr{Cvoid},) prob
    prob = CbcModel(p)
    finalizer(deleteModel, prob)
    return prob
end

function setParameter(prob::CbcModel, name::String, value::String)
    @assert isascii(name)
    @assert isascii(value)
    @cbc_ccall setParameter Cvoid (Ptr{Cvoid},Ptr{UInt8},Ptr{UInt8}) prob name value
end

# TODO: registerCallBack clearCallBack

function solve(prob::CbcModel)
    @cbc_ccall solve Cint (Ptr{Cvoid},) prob
end

@getproperty Float64 sumPrimalInfeasibilities
@getproperty Cint numberPrimalInfeasibilities

# TODO: checkSolution

@getproperty Cint getIterationCount

for s in (:isAbandoned, :isProvenOptimal, :isProvenInfeasible, :isContinuousUnbounded,
    :isNodeLimitReached, :isSecondsLimitReached, :isSolutionLimitReached, :isInitialSolveAbandoned, :isInitialSolveProvenOptimal, :isInitialSolveProvenPrimalInfeasible)
    @eval function ($s)(prob::CbcModel)
        check_problem(prob)
        v = @cbc_ccall $s Cint (Ptr{Cvoid},) prob
        return v != 0
    end
end

@getproperty Float64 getObjValue
@getproperty Float64 getBestPossibleObjValue
@getproperty Cint getNodeCount

# TODO: printSolution

function addSOS(prob::CbcModel, numRows::Integer, rowStarts::Vector{Cint},
    colIndices::Vector{Cint}, weights::Vector{Float64}, typ::Integer)

    @cbc_ccall(addSOS,Cvoid,(Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint},
                            Ptr{Float64}, Cint), prob, numRows,
                            rowStarts .- convert(Cint,1),
                            colIndices .- convert(Cint,1),
                            weights, typ)
end

# see Cbc_C_Interface.h documentation
@getproperty Cint status
@getproperty Cint secondaryStatus

# Both values are set to `-1` in the constructor and `resetModel`
function optimizeNotCalled(prob::CbcModel)
    return status(prob) == -1 && secondaryStatus(prob) == -1
end

end
