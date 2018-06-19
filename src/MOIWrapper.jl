export CbcOptimizer


using MathOptInterface
using Cbc.CbcCInterface

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const CbcCI = CbcCInterface

mutable struct CbcOptimizer <: MOI.AbstractOptimizer
    inner::CbcModel
    # env
    objFunctConstant::Float64
    CbcOptimizer() = new(CbcModel(), 0.0) # Initializes with an empty model
end

mutable struct CbcModelFormat
    nbRows::Int
    nbCols::Int
    constraint_matrix::SparseMatrixCSC{Float64,Int32}
    col_lb::Vector{Float64}
    col_ub::Vector{Float64}
    obj::Vector{Float64}
    row_lb::Vector{Float64}
    row_ub::Vector{Float64}
    function CbcModelFormat(nbRows::Int, nbCols::Int)
        obj = fill(0.0, nbCols)
        col_lb = fill(-Inf, nbCols)
        col_ub = fill(Inf, nbCols)
        row_lb = fill(-Inf, nbRows)
        row_ub = fill(Inf, nbRows)
        constraint_matrix = sparse(Int[], Int[], Float64[], nbRows, nbCols)
        new(nbRows, nbCols, constraint_matrix, col_lb, col_ub, obj, row_lb, row_ub)
    end
end


function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.SingleVariable, s::MOI.EqualTo)
    cbcModelFormat.col_lb[mapping.varmap[f.variable].value] = s.value
    cbcModelFormat.col_ub[mapping.varmap[f.variable].value] = s.value
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.SingleVariable, s::MOI.LessThan)
    cbcModelFormat.col_ub[mapping.varmap[f.variable].value] = s.upper
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.SingleVariable, s::MOI.GreaterThan)
    cbcModelFormat.col_lb[mapping.varmap[f.variable].value] = s.lower
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.SingleVariable, s::MOI.Interval)
    cbcModelFormat.col_lb[mapping.varmap[f.variable].value] = s.lower
    cbcModelFormat.col_ub[mapping.varmap[f.variable].value] = s.upper
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.ScalarAffineFunction, s::MOI.EqualTo)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.conmap[ci].value,mapping.varmap[term.variable_index].value] += term.coefficient
    end
    cbcModelFormat.row_lb[mapping.conmap[ci].value] = s.value - f.constant
    cbcModelFormat.row_ub[mapping.conmap[ci].value] = s.value - f.constant
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.ScalarAffineFunction, s::MOI.GreaterThan)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.conmap[ci].value,mapping.varmap[term.variable_index].value] += term.coefficient
    end
    cbcModelFormat.row_lb[mapping.conmap[ci].value] = s.lower - f.constant
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.ScalarAffineFunction, s::MOI.LessThan)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.conmap[ci].value, mapping.varmap[term.variable_index].value] += term.coefficient
    end
    cbcModelFormat.row_ub[mapping.conmap[ci].value] = s.upper - f.constant
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.ScalarAffineFunction, s::MOI.Interval)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.conmap[ci].value, mapping.varmap[term.variable_index].value] += term.coefficient
    end
    cbcModelFormat.row_ub[mapping.conmap[ci].value] = s.upper - f.constant
    cbcModelFormat.row_lb[mapping.conmap[ci].value] = s.lower - f.constant
end


function loadObj(cbcModelFormat::CbcModelFormat, mapping::MOIU.IndexMap,
    f::MOI.ScalarAffineFunction)
    for term in f.terms
        cbcModelFormat.obj[mapping.varmap[term.variable_index].value] += term.coefficient
    end
end

function copyconstraints(cbcModelFormat::CbcModelFormat, userOptimizer::MOI.ModelLike, mapping::MOIU.IndexMap, ::Type{F}, ::Type{S}) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    for ci in MOI.get(userOptimizer, MOI.ListOfConstraintIndices{F,S}())
        f = MOI.get(userOptimizer, MOI.ConstraintFunction(), ci)
        s = MOI.get(userOptimizer,  MOI.ConstraintSet(), ci)
        loadConstraint(ci, cbcModelFormat, mapping, f, s)
    end
end

```
    Create inner model cbcOptimizer based on abstract model userOptimizer provided by user.
    Fill the object of type CbcModelFormat:
        constraint_matrix::AbstractMatrix,
        col_lb::VecOrNothing,
        col_ub::VecOrNothing,
        obj::VecOrNothing,
        row_lb::VecOrNothing,
        row_ub::VecOrNothing,
    These are needed by function loadProblem of CbcCInterface.

```
function MOI.copy!(cbcOptimizer::CbcOptimizer,
    userOptimizer::MOI.ModelLike; copynames=false)

    mapping = MOIU.IndexMap()

    nbCols = MOI.get(userOptimizer, MOI.NumberOfVariables())
    varIndex = MOI.get(userOptimizer, MOI.ListOfVariableIndices())
    for i in 1:nbCols
        mapping.varmap[varIndex[i]] = MOI.VariableIndex(i)
    end

    zeroOneIndices = Int[]
    integerIndices = Int[]
    listOfConstraints = MOI.get(userOptimizer, MOI.ListOfConstraints())
    nbRows = 0
    for (F,S) in listOfConstraints
        if !(MOI.supportsconstraint(cbcOptimizer, F, S))
            return MOI.CopyResult(MOI.CopyUnsupportedConstraint, "Cbc MOI Interface does not support constraints of type " * (F,S) * ".", nothing)
        end

        ci = MOI.get(userOptimizer, MOI.ListOfConstraintIndices{F,S}())

        if S <: Union{MOI.ZeroOne, MOI.Integer}
            for i in 1:length(ci)
                f = MOI.get(userOptimizer, MOI.ConstraintFunction(), ci[i])
                S == MOI.ZeroOne && push!(zeroOneIndices, mapping.varmap[f.variable].value)
                S == MOI.Integer && push!(integerIndices, mapping.varmap[f.variable].value)
            end
        end

        ## Update conmap for (F,S) for F != MOI.SingleVariable
        ## Single variables are treated by bounds in Cbc, so no
        ## need to add a row
        if F != MOI.SingleVariable
            for i in 1:length(ci)
                mapping.conmap[ci[i]] = MOI.ConstraintIndex{F,S}(nbRows + i)
            end
            nbRows += MOI.get(userOptimizer, MOI.NumberOfConstraints{F,S}())
        end
    end

    cbcModelFormat = CbcModelFormat(nbRows, nbCols)

    for (F,S) in listOfConstraints
        if !(S <: Union{MOI.ZeroOne, MOI.Integer})
            copyconstraints(cbcModelFormat, userOptimizer, mapping, F, S)
        end
    end

    for colIdx in zeroOneIndices
        if cbcModelFormat.col_lb[colIdx] < 0.0
            cbcModelFormat.col_lb[colIdx] = 0.0
        end
        if cbcModelFormat.col_ub[colIdx] > 1.0
            cbcModelFormat.col_ub[colIdx] = 1.0
        end
    end


    ## Copy objective function
    objF = MOI.get(userOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    loadObj(cbcModelFormat, mapping, objF)
    cbcOptimizer.objFunctConstant = objF.constant
    sense = MOI.get(userOptimizer, MOI.ObjectiveSense())
    MOI.set!(cbcOptimizer, MOI.ObjectiveSense(), sense)

    ## Load the problem to Cbc
    CbcCI.loadProblem(cbcOptimizer.inner, cbcModelFormat.constraint_matrix, cbcModelFormat.col_lb, cbcModelFormat.col_ub, cbcModelFormat.obj, cbcModelFormat.row_lb, cbcModelFormat.row_ub)

    ## Set integer variables
    for colIdx in vcat(integerIndices, zeroOneIndices)
        CbcCI.setInteger(cbcOptimizer.inner, colIdx-1)
    end

    return MOI.CopyResult(MOI.CopySuccess, "Model was copied succefully.", MOIU.IndexMap(mapping.varmap, mapping.conmap))
end


function MOI.optimize!(cbcOptimizer::CbcOptimizer)
    # Call solve function
    CbcCI.solve(cbcOptimizer.inner)
end



## canadd, canset, canget functions

function MOI.canaddvariable(cbcOptimizer::CbcOptimizer)
    return false
end

## supports constraints


MOI.supportsconstraint(::CbcOptimizer, ::Type{<:Union{MOI.ScalarAffineFunction{Float64}, MOI.SingleVariable}}, ::Type{<:Union{MOI.EqualTo{Float64}, MOI.Interval{Float64}, MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.ZeroOne, MOI.Integer}}) = true


MOI.supports(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true



## Set functions


# empty!
function MOI.empty!(cbcOptimizer::CbcOptimizer)
    cbcOptimizer.inner = CbcModel()
end


function MOI.set!(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    sense == MOI.MinSense && CbcCI.setObjSense(cbcOptimizer.inner, 1)
    sense == MOI.MaxSense && CbcCI.setObjSense(cbcOptimizer.inner, -1)
end


## Get functions


function MOI.canget(cbcOptimizer::CbcOptimizer, object::MOI.PrimalStatus)
    object.N != 1 && return false
    return MOI.get(cbcOptimizer, MOI.ResultCount()) == 1
end


MOI.canget(cbcOptimizer::CbcOptimizer, object::Union{MOI.NodeCount, MOI.ResultCount, MOI.TerminationStatus, MOI.ObjectiveSense, MOI.ObjectiveValue, MOI.ObjectiveBound, MOI.NumberOfVariables}) = true

MOI.canget(cbcOptimizer::CbcOptimizer, object::MOI.VariablePrimal, indexTypeOrObject::Type{MOI.VariableIndex}) = true


MOI.isempty(cbcOptimizer::CbcOptimizer) = (CbcCI.getNumCols(cbcOptimizer.inner) == 0 && CbcCI.getNumRows(cbcOptimizer.inner) == 0)

MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.NumberOfVariables) = getNumCols(cbcOptimizer.inner)

MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveValue) = (CbcCI.getObjValue(cbcOptimizer.inner) + cbcOptimizer.objFunctConstant)

MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveBound) = CbcCI.getBestPossibleObjValue(cbcOptimizer.inner)

MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.NodeCount) = CbcCI.getNodeCount(cbcOptimizer.inner)

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.VariablePrimal, ref::MOI.VariableIndex)
    variablePrimals = CbcCI.getColSolution(cbcOptimizer.inner)
    return variablePrimals[ref.value]
end

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.VariablePrimal, ref::Vector{MOI.VariableIndex})
    variablePrimals = CbcCI.getColSolution(cbcOptimizer.inner)
    return [variablePrimals[vi.value] for vi in ref]
end


function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.ResultCount)
    isProvenInfeasible(cbcOptimizer.inner) && return 0
    isContinuousUnbounded(cbcOptimizer.inner) && return 0
    isAbandoned(cbcOptimizer.inner) && return 0
    CbcCI.getObjValue(cbcOptimizer.inner) >= 1e300 && return 0
    return 1
end


function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveSense)
    CbcCI.getObjSense(cbcOptimizer.inner) == 1 && return MOI.MinSense
    CbcCI.getObjSense(cbcOptimizer.inner) == -1 && return MOI.MaxSense
end



function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.TerminationStatus)

    isProvenInfeasible(cbcOptimizer.inner) && return MOI.InfeasibleNoResult
    isContinuousUnbounded(cbcOptimizer.inner) && return MOI.InfeasibleOrUnbounded
    isNodeLimitReached(cbcOptimizer.inner) && return MOI.NodeLimit
    isSecondsLimitReached(cbcOptimizer.inner) && return MOI.TimeLimit
    isSolutionLimitReached(cbcOptimizer.inner) && return MOI.SolutionLimit
    isProvenOptimal(cbcOptimizer.inner) && return MOI.Success
    isInitialSolveProvenOptimal(cbcOptimizer.inner) && return MOI.Success
    MOI.get(cbcOptimizer, MOI.ResultCount()) == 1 && return MOI.Success
    isAbandoned(cbcOptimizer.inner) && return MOI.Interrupted

    error("Internal error: Unrecognized solution status")
end

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.PrimalStatus)
    isProvenOptimal(cbcOptimizer.inner) && return MOI.FeasiblePoint
    isInitialSolveProvenOptimal(cbcOptimizer.inner) && return MOI.FeasiblePoint
    isProvenInfeasible(cbcOptimizer.inner) && return MOI.InfeasiblePoint
end




#
