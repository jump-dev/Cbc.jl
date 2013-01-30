

module CoinMP

include(joinpath(Pkg.dir(),"CoinMP","deps","ext.jl"))


const CONTINUOUS = 1
const INTEGER = 2 # including binary

export CONTINUOUS, INTEGER,
    mixintprog,
    CoinProblem,
    LoadMatrix,
    LoadInteger,
    GetColCount,
    GetRowCount,
    OptimizeProblem,
    GetSolutionStatus,
    GetSolutionText,
    GetSolutionValues,
    GetObjectValue

## High-level interface
# (z, x, flag) = mixintprog(f, A, rowlb, rowub, lb, ub, vartype)
# Solves: 
#   z = min_{x} (f' * x)
#
# where the vector x is subject to these constraints:
#
#   rowlb <= A * x <= rowub
#   lb <= x <= ub
#
# The return flag is a string describing the status of the optimization.
# Success is "Optimal solution found"
# The following arguments may be nothing or [], in which case 
# the default values are vectors of:
# rowlb -- -Inf
# rowub -- Inf
# collb -- 0
# colub -- Inf
#
# vartype is a vector specifying the type of each variable
# These can be CONTINUOUS or INTEGER

function mixintprog(f, A, rowlb, rowub, lb, ub, vartype::Vector)
    c = CoinProblem()
    LoadMatrix(c, 1, 0., f,  lb, ub, rowlb, rowub, A)
    ncol = GetColCount(c)
    if length(vartype) != ncol
        error("Length of vector must equal number of variables in the problem")
    end
    coltype = Array(Uint8,ncol)
    for i in 1:ncol
        if vartype[i] == INTEGER
            coltype[i] = 'I'
        else
            coltype[i] = 'C'
        end
    end
    LoadInteger(c, coltype)
    OptimizeProblem(c)
    stat = GetSolutionText(c)
    if stat != "Optimal solution found"
        return (nothing, nothing, stat)
    else
        return (GetObjectValue(c),GetSolutionValues(c),stat)
    end
end

# helper macros/functions

macro coin_ccall(func, args...)
    f = "Coin$(func)"
    quote
        ccall(($f,coinmp_lib), $(args...))
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
    @coin_ccall UnloadProblem Int32 (Ptr{Void},) prob.p
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
    if isequal(a, nothing) || isa(a, Array{None})
        return C_NULL
    else # todo: helpful message if convert fails
        if length(a) != len
            error("Expected vector to have length $len")
        end
        return convert(Vector{T},a)
    end
end

# actual interface

#SOLVAPI int SOLVCALL CoinLoadMatrix(HPROB hProb, 
#				int ColCount, int RowCount, int NZCount, int RangeCount, 
#				int ObjectSense, double ObjectConst, double* ObjectCoeffs, 
#				double* LowerBounds, double* UpperBounds, const char* RowType, 
#				double* RHSValues, double* RangeValues, int* MatrixBegin, 
#				int* MatrixCount, int* MatrixIndex, double* MatrixValues);


function LoadMatrix(prob::CoinProblem, objective_sense::Integer,
    objective_offset::Float64, objective_coeffs::VecOrNothing, 
    col_lb::VecOrNothing, col_ub::VecOrNothing,
    row_lb::VecOrNothing, row_ub::VecOrNothing,
    constraint_matrix::AbstractMatrix{Float64})
    check_problem(prob)

    mat = convert(SparseMatrixCSC{Float64,Int32},convert(SparseMatrixCSC,constraint_matrix))
    nrow,ncol = mat.m, mat.n

    ret = @coin_ccall LoadMatrix Int32 (Ptr{Void}, Int32, Int32, Int32, Int32, Int32,
        Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Uint8}, Ptr{Float64},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, 
        Ptr{Float64}) prob.p ncol nrow nnz(mat) 0 objective_sense objective_offset vec_or_null(Float64, objective_coeffs, ncol) vec_or_null(Float64, col_lb, ncol) vec_or_null(Float64, col_ub, ncol) C_NULL vec_or_null(Float64, row_lb, nrow) vec_or_null(Float64, row_ub, nrow) mat.colptr-int32(1) C_NULL mat.rowval-int32(1) mat.nzval
    if ret != 0
        error("Internal error in CoinLoadMatrix")
    end
end

# 'B' - binary
# 'I' - integer
# otherwise continuous
function LoadInteger(prob::CoinProblem, column_type::Vector{Uint8})
    check_problem(prob)
    if length(column_type) != GetColCount(prob)
        error("Length of vector must equal number of variables in the problem")
    end
    ret = @coin_ccall LoadInteger Int32 (Ptr{Void}, Ptr{Uint8}) prob.p column_type
    if ret != 0
        error("Internal error in CoinLoadInteger")
    end
end

   
function GetColCount(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetColCount Int32 (Ptr{Void},) prob.p
end

function GetRowCount(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetRowCount Int32 (Ptr{Void},) prob.p
end

function OptimizeProblem(prob::CoinProblem)
    check_problem(prob)
    # Method argument is ignored?
    ret = @coin_ccall OptimizeProblem Int32 (Ptr{Void},Int32) prob.p 0
    if ret != 0
        error("Internal error in CoinOptimizeProblem")
    end
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

# return vector of primal solutions
# we could return more in the future
function GetSolutionValues(prob::CoinProblem)
    check_problem(prob)
    ncol = GetColCount(prob)
    x = Array(Float64,ncol)
    ret = @coin_ccall GetSolutionValues Int32 (Ptr{Void},
        Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}) prob.p x C_NULL C_NULL C_NULL
    if ret != 0
        error("Internal error in CoinGetSolutionValues")
    end
    return x
end

function GetObjectValue(prob::CoinProblem)
    check_problem(prob)
    @coin_ccall GetObjectValue Float64 (Ptr{Void},) prob.p
end


end # module
