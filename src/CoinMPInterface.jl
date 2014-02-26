module CoinMPInterface

import Cbc

const CONTINUOUS = 1
const INTEGER = 2 # including binary

export CONTINUOUS, INTEGER,
    optionList,
    setOption,
    getOptionValue,
    getOptionDefault,
    getOptionMin,
    getOptionMax,
    CoinProblem,
    GetSolverName,
    GetSolverVersionStr,
    GetSolverVersion,
    GetInfinity,
    LoadMatrix,
    LoadInteger,
    LoadPriority,
    LoadSemiCont,
    CheckProblem,
    GetColCount,
    GetRowCount,
    RegisterMsgLogCallback,
    RegisterLPIterCallback,
    RegisterMipNodeCallback,
    OptimizeProblem,
    GetSolutionStatus,
    GetSolutionText,
    GetObjectValue,
    GetMipBestBound,
    GetIterCount,
    GetMipNodeCount,
    GetSolutionValues,
    ReadFile,
    WriteFile,
    OpenLogFile,
    CloseLogFile

# helper macros/functions

macro coin_ccall(func, args...)
    f = "Coin$(func)"
    quote
        @unix_only ret = ccall(($f,Cbc.libcoinmp), $(args...))
        @windows_only ret = ccall(($f,Cbc.libcoinmp), stdcall, $(args...))
        ret
    end
end

macro coin_checkedccall(func, args...)
    f = "Coin$(func)"
    quote
        @unix_only ret = ccall(($f,Cbc.libcoinmp), Int32, $(args...))
        @windows_only ret = ccall(($f,Cbc.libcoinmp), stdcall, Int32, $(args...))
        if ret != 0
            error("Internal error in ", $f)
        end
    end
end


type CoinProblem
    p::Ptr{Void}
    function CoinProblem()
        p = @coin_ccall CreateProblem Ptr{Void} (Ptr{Uint8},) bytestring("")
        prob = new(p)
        finalizer(prob, delete_problem)
        return prob
    end
end

function delete_problem(prob::CoinProblem)
    if prob.p == C_NULL
        return
    end
    @coin_checkedccall UnloadProblem (Ptr{Void},) prob.p
    prob.p = C_NULL
    return
end

function check_problem(prob::CoinProblem)
   if prob.p == C_NULL
       error("Invalid CoinProblem")
   end
   return true
end

# copied from Clp.jl
typealias VecOrNothing Union(Vector,Nothing)
function vec_or_null{T}(::Type{T}, a::VecOrNothing, len::Integer)
    if isequal(a, nothing) || isa(a, Array{None}) || length(a) == 0
        return C_NULL
    else # todo: helpful message if convert fails
        if length(a) != len
            error("Expected vector to have length $len")
        end
        return convert(Vector{T},a)
    end
end

## Convenience methods for options
# There's no documentation for the meanings of various options and their values
# The only way is to dig into the CoinMP source or find examples...

# returns the list of options
function optionList()
    dummy = CoinProblem()
    l = ASCIIString[]
    for i in 1:200
        s = GetOptionName(dummy,i)
        if s != ""
            push!(l,s)
        end
    end
    return l
end

function setOption(prob::CoinProblem, name::ASCIIString, value)
    id = LocateOptionName(prob,name)
    t = GetOptionType(prob,id)
    if t == 2 || t == 1
        SetIntOption(prob, id, value)
    elseif t == 4
        SetRealOption(prob, id, value)
    else
        error("Unkown option type, this is a bug")
    end
end

for (n1, n2) in ( ("Value",""), ("Default","DefaultValue"),
    ("Min","MinValue"), ("Max","MaxValue") )
    fname = symbol(string("getOption",n1))
    coinint = symbol(string("GetIntOption",n2))
    coinreal = symbol(string("GetRealOption",n2))
    @eval begin function $(fname)(prob::CoinProblem, name::ASCIIString)
        id = LocateOptionName(prob,name)
        t = GetOptionType(prob,id)
        if t == 2 || t == 1
            $coinint(prob, id)
        elseif t == 4
            $coinreal(prob, id)
        else
            error("Unkown option type, this is a bug")
        end
    end end
end


## actual interface

# CoinInitSolver and CoinFreeSolver don't do anything

function GetSolverName()
    s = @coin_ccall GetSolverName Ptr{Uint8} ()
    bytestring(s)
end

function GetSolverVersionStr()
    s = @coin_ccall GetSolverVersionStr Ptr{Uint8} ()
    bytestring(s)
end

function GetSolverVersion()
    @coin_ccall GetSolverVersion Float64 ()
end

# not sure the point of CoinGetFeatures and CoinGetMethods

function GetInfinity()
    @coin_ccall GetInfinity Float64 ()
end

function LoadMatrix(prob::CoinProblem, objective_sense::Integer,
    objective_offset::Float64, objective_coeffs::VecOrNothing,
    col_lb::VecOrNothing, col_ub::VecOrNothing,
    row_lb::VecOrNothing, row_ub::VecOrNothing,
    constraint_matrix::AbstractMatrix)
    check_problem(prob)

    mat = convert(SparseMatrixCSC{Float64,Int32},constraint_matrix)
    nrow,ncol = mat.m, mat.n

    @coin_checkedccall LoadMatrix (Ptr{Void}, Int32, Int32, Int32, Int32, Int32,
        Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Uint8}, Ptr{Float64},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}) prob.p ncol nrow nnz(mat) 0 objective_sense objective_offset vec_or_null(Float64, objective_coeffs, ncol) vec_or_null(Float64, col_lb, ncol) vec_or_null(Float64, col_ub, ncol) C_NULL vec_or_null(Float64, row_lb, nrow) vec_or_null(Float64, row_ub, nrow) mat.colptr-int32(1) C_NULL vec_or_null(Int32, mat.rowval-int32(1),nnz(mat)) vec_or_null(Float64, mat.nzval, nnz(mat))
end

# TODO: CoinLoadNames

# CoinLoadInitValues doesn't do anything

# 'B' - binary
# 'I' - integer
# 'C' - continuous
function LoadInteger(prob::CoinProblem, column_type::Vector{Uint8})
    check_problem(prob)
    if length(column_type) != GetColCount(prob)
        error("Length of vector must equal number of variables in the problem")
    end
    @coin_checkedccall LoadInteger (Ptr{Void}, Ptr{Uint8}) prob.p column_type
end

function LoadPriority(prob::CoinProblem,
    PriorIndex::Vector{Int32}, PriorValues::Vector{Int32})
    # PriorBranch argument is ignored?
    check_problem(prob)
    len = length(PriorIndex)
    PriorIndex = PriorIndex - int32(1) # input is 1-based
    if len != length(PriorValues)
        error("All input vectors must be same length")
    end
    if len > GetColCount(prob)
        error("More priorities than variables")
    end
    @coin_checkedccall LoadPriority (Ptr{Void}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}) prob.p PriorIndex PriorValues C_NULL
end

# TODO: CoinLoadSos

function LoadSemiCont(prob::CoinProblem, SemiIndex::Vector{Int32})
    check_problem(prob)
    SemiIndex = SemiIndex - int32(1)
    @coin_checkedccall LoadSemiCont (Ptr{Void}, Int32, Ptr{Int32}) prob.p length(SemiIndex) SemiIndex
end

# CoinLoadQuadratic and CoinLoadNonlinear are not implemented by CoinMP

# This returns 0 if the problem data is ok, otherwise it returns
# one of the SOLV_CHECK_* error codes defined in CoinMP.h
function CheckProblem(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall CheckProblem Int32 (Ptr{Void},) prob.p
end


# TODO: CoinGetProblemName -- is this useful?

function GetColCount(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetColCount Int32 (Ptr{Void},) prob.p
end

function GetRowCount(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetRowCount Int32 (Ptr{Void},) prob.p
end

# TODO: CoinGetColName, CoinGetRowName

# Pass in c-pointer to callback function. (from cfunction())
# See CoinMP.h for signature
# TODO: test these
function RegisterMsgLogCallback(prob::CoinProblem, func::Ptr{Void})
    check_problem(prob)
    @coin_checkedccall RegisterMsgLogCallback (Ptr{Void},Ptr{Void},Ptr{Void}) prob.p func C_NULL
end

function RegisterLPIterCallback(prob::CoinProblem, func::Ptr{Void})
    check_problem(prob)
    @coin_checkedccall RegisterLPIterCallback (Ptr{Void},Ptr{Void},Ptr{Void}) prob.p func C_NULL
end

function RegisterMipNodeCallback(prob::CoinProblem, func::Ptr{Void})
    check_problem(prob)
    @coin_checkedccall RegisterMipNodeCallback (Ptr{Void},Ptr{Void},Ptr{Void}) prob.p func C_NULL
end

function OptimizeProblem(prob::CoinProblem)
    check_problem(prob)
    # Method argument is ignored?
    @coin_checkedccall OptimizeProblem (Ptr{Void},Int32) prob.p 0
end

# no documentation for what these values mean...
function GetSolutionStatus(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetSolutionStatus Int32 (Ptr{Void},) prob.p
end

function GetSolutionText(prob::CoinProblem)
    check_problem(prob)
    s = @coin_ccall GetSolutionText Ptr{Uint8} (Ptr{Void},) prob.p
    return bytestring(s)
end

function GetObjectValue(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetObjectValue Float64 (Ptr{Void},) prob.p
end

function GetMipBestBound(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetMipBestBound Float64 (Ptr{Void},) prob.p
end

function GetIterCount(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetIterCount Int32 (Ptr{Void},) prob.p
end

function GetMipNodeCount(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetMipNpdeCount Int32 (Ptr{Void},) prob.p
end

# return vector of primal solutions
# we could return more in the future
function GetSolutionValues(prob::CoinProblem)
    check_problem(prob)
    ncol = GetColCount(prob)
    x = Array(Float64,ncol)
    @coin_checkedccall GetSolutionValues (Ptr{Void},
        Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}) prob.p x C_NULL C_NULL C_NULL
    return x
end

# TODO: GetSolutionRanges -- what does this mean? sensitivity analysis?

# TODO: GetSolutionBasis -- pretty useless since you can't input a starting basis

# the only file type currently supported is MPS (3)
function ReadFile(prob::CoinProblem, FileType::Integer, FileName::ASCIIString)
    check_problem(prob)
    @coin_checkedccall ReadFile (Ptr{Void},Int32,Ptr{Uint8}) prob.p FileType FileName
end

function WriteFile(prob::CoinProblem, FileType::Integer, FileName::ASCIIString)
    check_problem(prob)
    @coin_checkedccall WriteFile (Ptr{Void},Int32,Ptr{Uint8}) prob.p FileType FileName
end

function OpenLogFile(prob::CoinProblem, FileName::ASCIIString)
    check_problem(prob)
    @coin_checkedccall OpenLogFile (Ptr{Void},Ptr{Uint8}) prob.p FileName
end

function CloseLogFile(prob::CoinProblem)
    check_problem(prob)
    @coin_checkedccall CloseLogFile (Ptr{Void},) prob.p
end

function GetOptionCount(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetOptionCount Int32 (Ptr{Void},) prob.p
end

function LocateOptionID(prob::CoinProblem, OptionID::Integer)
    check_problem(prob)
    @coin_ccall LocateOptionID Int32 (Ptr{Void},Int32) prob.p OptionID
end

function LocateOptionName(prob::CoinProblem, OptionName::ASCIIString)
    check_problem(prob)
    @coin_ccall LocateOptionName Int32 (Ptr{Void},Ptr{Uint8}) prob.p OptionName
end

# Don't know what this does
function GetOptionID(prob::CoinProblem, OptionNr::Integer)
    check_problem(prob)
    @coin_ccall GetOptionID Int32 (Ptr{Void},Int32) prob.p OptionNr
end

# TODO: CoinGetOptionInfo, CoinGetIntOptionMinMax, CoinGetRealOptionMinMax

# TODO: CoinGetOptionGroup

# 1 -- 0/1, 2 -- INT, 4 -- REAL,
function GetOptionType(prob::CoinProblem, OptionID::Integer)
    check_problem(prob)
    @coin_ccall GetOptionType Int32 (Ptr{Void},Int32) prob.p OptionID
end

for s in (:GetIntOptionDefaultValue, :GetIntOptionMinValue, :GetIntOptionMaxValue)
    @eval begin
        function $(s)(prob::CoinProblem, OptionID::Integer)
            check_problem(prob)
            @coin_ccall $s Int32 (Ptr{Void},Int32) prob.p OptionID
        end
    end
end

for s in (:GetRealOptionDefaultValue, :GetRealOptionMinValue, :GetRealOptionMaxValue)
    @eval begin
        function $(s)(prob::CoinProblem, OptionID::Real)
            check_problem(prob)
            @coin_ccall $s Float64 (Ptr{Void},Int32) prob.p OptionID
        end
    end
end

for s in (:GetOptionName, :GetOptionShortName, :GetStringOption)
    @eval begin
        function $(s)(prob::CoinProblem, OptionID::Integer)
            check_problem(prob)
            bytestring(@coin_ccall GetOptionName Ptr{Uint8} (Ptr{Void},Int32) prob.p OptionID)
        end
    end
end

# TODO: CoinGetOptionChanged

function GetIntOption(prob::CoinProblem, OptionID::Integer)
    check_problem(prob)
    @coin_ccall GetIntOption Int32 (Ptr{Void},Int32) prob.p OptionID
end

function SetIntOption(prob::CoinProblem, OptionID::Integer, Value::Integer)
    check_problem(prob)
    @coin_checkedccall SetIntOption (Ptr{Void},Int32,Int32) prob.p OptionID Value
end

function GetRealOption(prob::CoinProblem, OptionID::Integer)
    check_problem(prob)
    @coin_ccall GetRealOption Float64 (Ptr{Void},Int32) prob.p OptionID
end

function SetRealOption(prob::CoinProblem, OptionID::Integer, Value::Real)
    check_problem(prob)
    @coin_checkedccall SetRealOption (Ptr{Void},Int32,Float64) prob.p OptionID Value
end

function SetStringOption(prob::CoinProblem, OptionID::Integer, Value::ASCIIString)
    check_problem(prob)
    @coin_checkedccall SetStringOption (Ptr{Void},Int32,Ptr{Uint8}) prob.p OptionID Value
end

end
