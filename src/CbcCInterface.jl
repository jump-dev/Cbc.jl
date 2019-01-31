module CbcCInterface
using Compat: Cvoid, @compat
using Compat.SparseArrays

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
    args = map(esc,args)
    f = "Cbc_$(func)"
    quote
        ccall(($f,Cbc.libcbcsolver), $(args...))
    end
end

mutable struct CbcModel
    p::Ptr{Cvoid}
    function CbcModel()
        p = @cbc_ccall newModel Ptr{Cvoid} ()
        prob = new(p)
        @compat finalizer(deleteModel, prob)
        return prob
    end
end

function deleteModel(prob::CbcModel)
    if prob.p == C_NULL
        return
    end
    GC.@preserve prob begin
        @cbc_ccall deleteModel Cvoid (Ptr{Cvoid},) prob.p
    end
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
        GC.@preserve prob begin
            @cbc_ccall $name $T (Ptr{Cvoid},) prob.p
        end
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

    GC.@preserve prob mat begin
        @cbc_ccall loadProblem Cvoid (Ptr{Cvoid}, Int32, Int32, Ptr{CoinBigIndex},
            Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}) prob.p ncol nrow mat.colptr.-Int32(1) mat.rowval.-Int32(1) mat.nzval vec_or_null(Float64, col_lb, ncol) vec_or_null(Float64, col_ub, ncol) vec_or_null(Float64, obj, ncol) vec_or_null(Float64, row_lb, nrow) vec_or_null(Float64, row_ub, nrow)
    end
end

function readMps(prob::CbcModel, filename::String)
    check_problem(prob)
    @assert isascii(filename)
    GC.@preserve prob filename begin
        @cbc_ccall readMps Cint (Ptr{Cvoid}, Ptr{UInt8}) prob.p filename
    end
end

function writeMps(prob::CbcModel, filename::String)
    check_problem(prob)
    @assert isascii(filename)
    GC.@preserve prob filename begin
        @cbc_ccall writeMps Cint (Ptr{Cvoid}, Ptr{UInt8}) prob.p filename
    end
end

function setInitialSolution(prob::CbcModel, array::Vector{Float64})
    check_problem(prob)
    GC.@preserve prob array begin
        @cbc_ccall setInitialSolution Cvoid (Ptr{Cvoid}, Ptr{Float64}) prob.p array
    end
end

function problemName(prob::CbcModel)
    check_problem(prob)
    a = Array(UInt8, 100)
    GC.@preserve prob a begin
        @cbc_ccall problemName Cvoid (Ptr{Cvoid},Cint,Ptr{UInt8}) prob.p 100 a
    end
    return string(a)
end

#function setProblemName(prob::CbcModel)
#    check_problem(prob)
#    @cbc_ccall setProblemName Cint (Ptr{Cvoid},Ptr{UInt8}) prob.p bytestring(a)
#end

function getNumElements(prob::CbcModel)
    check_problem(prob)
    GC.@preserve prob begin
        @cbc_ccall getNumElements Cint (Ptr{Cvoid},) prob.p
    end
end

function getVectorStarts(prob::CbcModel)
    check_problem(prob)
    GC.@preserve prob begin
        p = @cbc_ccall getVectorStarts Ptr{CoinBigIndex} (Ptr{Cvoid},) prob.p
    end
    num_cols = Int(getNumCols(prob))
    return copy(unsafe_wrap(Array,p,(num_cols+1,)))
end

function getIndices(prob::CbcModel)
    check_problem(prob)
    GC.@preserve prob begin
        p = @cbc_ccall getIndices Ptr{Cint} (Ptr{Cvoid},) prob.p
    end
    nnz = Int(getNumElements(prob))
    return copy(unsafe_wrap(Array,p,(nnz,)))
end

function getElements(prob::CbcModel)
    check_problem(prob)
    GC.@preserve prob begin
        p = @cbc_ccall getElements Ptr{Float64} (Ptr{Cvoid},) prob.p
    end
    nnz = Int(getNumElements(prob))
    return copy(unsafe_wrap(Array,p,(nnz,)))
end

# TODO: maxNameLength getRowName getColName setColName setRowName

@getproperty Cint getNumRows
@getproperty Cint getNumCols


# 1 : minimize, -1 : maximize
function setObjSense(prob::CbcModel, sense)
    check_problem(prob)
    GC.@preserve prob sense begin
        @cbc_ccall setObjSense Cvoid (Ptr{Cvoid}, Float64) prob.p sense
    end
end

@getproperty Float64 getObjSense

for s in (:getRowLower, :getRowUpper, :getRowActivity)
    @eval function ($s)(prob::CbcModel)
        check_problem(prob)
        nrow = Int(getNumRows(prob))
        GC.@preserve prob begin
            p = @cbc_ccall $s Ptr{Float64} (Ptr{Cvoid},) prob.p
        end
        return copy(unsafe_wrap(Array,p,(nrow,)))
    end
end

for s in (:getColLower, :getColUpper, :getObjCoefficients, :getColSolution)
    @eval function ($s)(prob::CbcModel)
        check_problem(prob)
        ncol = Int(getNumCols(prob))
        GC.@preserve prob begin
            p = @cbc_ccall $s Ptr{Float64} (Ptr{Cvoid},) prob.p
        end
        return copy(unsafe_wrap(Array,p,(ncol,)))
    end
end

for s in (:setRowUpper, :setRowLower, :setObjCoeff, :setColLower, :setColUpper)
    @eval function($s)(prob::CbcModel, index::Integer, value::Float64)
        check_problem(prob)
        GC.@preserve prob index value begin
            @cbc_ccall $s Cvoid (Ptr{Cvoid}, Cint, Float64) prob.p index value
        end
    end
end

function isInteger(prob::CbcModel, index::Integer)
    check_problem(prob)
    GC.@preserve prob index begin
        v = @cbc_ccall isInteger Cint (Ptr{Cvoid},Cint) prob.p index
    end
    return v == 1
end

function setContinuous(prob::CbcModel, index::Integer)
    check_problem(prob)
    GC.@preserve prob index begin
        @cbc_ccall setContinuous Cvoid (Ptr{Cvoid},Cint) prob.p index
    end
end

function setInteger(prob::CbcModel, index::Integer)
    check_problem(prob)
    GC.@preserve prob index begin
        @cbc_ccall setInteger Cvoid (Ptr{Cvoid},Cint) prob.p index
    end
end

function Base.copy(prob::CbcModel)
    GC.@preserve prob begin
        p = @cbc_ccall clone Ptr{Cvoid} (Ptr{Cvoid},) prob.p
    end
    prob = CbcModel(p)
    @compat finalizer(deleteModel, prob)
    return prob
end

function setParameter(prob::CbcModel, name::String, value::String)
    @assert isascii(name)
    @assert isascii(value)
    GC.@preserve prob name value begin
        @cbc_ccall setParameter Cvoid (Ptr{Cvoid},Ptr{UInt8},Ptr{UInt8}) prob.p name value
    end
end

# TODO: registerCallBack clearCallBack

function solve(prob::CbcModel)
    GC.@preserve prob begin
        @cbc_ccall solve Cint (Ptr{Cvoid},) prob.p
    end
end

@getproperty Float64 sumPrimalInfeasibilities
@getproperty Cint numberPrimalInfeasibilities

# TODO: checkSolution

@getproperty Cint getIterationCount

for s in (:isAbandoned, :isProvenOptimal, :isProvenInfeasible, :isContinuousUnbounded,
    :isNodeLimitReached, :isSecondsLimitReached, :isSolutionLimitReached, :isInitialSolveAbandoned, :isInitialSolveProvenOptimal, :isInitialSolveProvenPrimalInfeasible)
    @eval function ($s)(prob::CbcModel)
        check_problem(prob)
        GC.@preserve prob begin
            v = @cbc_ccall $s Cint (Ptr{Cvoid},) prob.p
        end
        return v != 0
    end
end

@getproperty Float64 getObjValue
@getproperty Float64 getBestPossibleObjValue
@getproperty Cint getNodeCount

# TODO: printSolution

function addSOS(prob::CbcModel, numRows::Integer, rowStarts::Vector{Cint},
    colIndices::Vector{Cint}, weights::Vector{Float64}, typ::Integer)

    GC.@preserve prob numRows rowStarts colIndices weights typ begin
        @cbc_ccall(addSOS,Cvoid,(Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint},
                                Ptr{Float64}, Cint),prob.p,numRows,
                                rowStarts .- convert(Cint,1),
                                colIndices .- convert(Cint,1),
                                weights, typ)
    end
end

# see Cbc_C_Interface.h documentation
@getproperty Cint status
@getproperty Cint secondaryStatus

# Both values are set to `-1` in the constructor and `resetModel`
function optimizeNotCalled(prob::CbcModel)
    return status(prob) == -1 && secondaryStatus(prob) == -1
end

end
