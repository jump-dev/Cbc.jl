export CbcOptimizer, copy!, optimize!, get


using MathOptInterface
using Cbc.CbcCInterface

const MOI = MathOptInterface

const SUPPORTED_OBJECTIVES = [
    MOI.ScalarAffineFunction{Float64}
]

const SUPPORTED_CONSTRAINTS = [
    (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}),
    (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}),
    (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}),
    (MOI.SingleVariable, MOI.EqualTo{Float64}),
    (MOI.SingleVariable, MOI.LessThan{Float64}),
    (MOI.SingleVariable, MOI.GreaterThan{Float64}),
    (MOI.SingleVariable, MOI.Interval{Float64}),
    (MOI.SingleVariable, MOI.ZeroOne),
    (MOI.SingleVariable, MOI.Integer)
]

mutable struct VarAndConstrMapping
    variable_mapping::Dict{MOI.VariableIndex, Int}
    constraint_mapping::Dict{MOI.ConstraintIndex, Int}
    VarAndConstrMapping() = new(Dict{MOI.VariableIndex, Int}(), Dict{MOI.ConstraintIndex, Int}())
end

mutable struct CbcOptimizer <: MOI.AbstractOptimizer
    inner::CbcModel
    mapping::VarAndConstrMapping
    # env
    params::Dict{Symbol,Any}
    CbcOptimizer() = new(CbcModel(), VarAndConstrMapping()) # Initializes with an empty model
end

mutable struct CbcModelFormat
    nbRows::Int
    nbCols::Int
    constraint_matrix::AbstractMatrix
    col_lb::Vector{Float64}
    col_ub::Vector{Float64}
    obj::Vector{Float64}
    row_lb::Vector{Float64}
    row_ub::Vector{Float64}
    function CbcModelFormat(nbRows::Int, nbCols::Int)
        col_lb = Vector{Float64}(nbCols)
        col_ub = Vector{Float64}(nbCols)
        obj = Vector{Float64}(nbCols)
        row_lb = Vector{Float64}(nbRows)
        row_ub = Vector{Float64}(nbRows)
        constraint_matrix = Array{Float64, 2}(nbRows, nbCols)
        for rowIdx in 1:nbRows
            row_lb[rowIdx] = -Inf
            row_ub[rowIdx] = Inf
            for colIdx in 1:nbCols
                constraint_matrix[rowIdx, colIdx] = 0.0
            end
        end
        for colIdx in 1:nbCols
            col_lb[colIdx] = -Inf
            col_ub[colIdx] = Inf
            obj[colIdx] = 0.0
        end

        new(nbRows, nbCols, constraint_matrix, col_lb, col_ub, obj, row_lb, row_ub)
    end
end


function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.SingleVariable, s::MOI.EqualTo)
    cbcModelFormat.col_lb[mapping.variable_mapping[f.variable]] = s.value
    cbcModelFormat.col_ub[mapping.variable_mapping[f.variable]] = s.value
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.SingleVariable, s::MOI.LessThan)
    cbcModelFormat.col_ub[mapping.variable_mapping[f.variable]] = s.upper
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.SingleVariable, s::MOI.GreaterThan)
    cbcModelFormat.col_lb[mapping.variable_mapping[f.variable]] = s.lower
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.SingleVariable, s::MOI.Interval)
    cbcModelFormat.col_lb[mapping.variable_mapping[f.variable]] = s.lower
    cbcModelFormat.col_ub[mapping.variable_mapping[f.variable]] = s.upper
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.SingleVariable, s::MOI.ZeroOne)
    cbcModelFormat.col_lb[mapping.variable_mapping[f.variable]] = 0
    cbcModelFormat.col_ub[mapping.variable_mapping[f.variable]] = 1
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.ScalarAffineFunction, s::MOI.EqualTo)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.constraint_mapping[ci]][mapping.variable_mapping[term.variable_index]] = term.coefficient
    end
    cbcModelFormat.row_lb[mapping.constraint_mapping[ci]] = s.value - f.constant
    cbcModelFormat.row_ub[mapping.constraint_mapping[ci]] = s.value - f.constant
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.ScalarAffineFunction, s::MOI.GreaterThan)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.constraint_mapping[ci]][mapping.variable_mapping[term.variable_index]] = term.coefficient
    end
    cbcModelFormat.row_lb[mapping.constraint_mapping[ci]] = s.lower - f.constant
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.ScalarAffineFunction, s::MOI.LessThan)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.constraint_mapping[ci]][mapping.variable_mapping[term.variable_index]] = term.coefficient
    end
    cbcModelFormat.row_ub[mapping.constraint_mapping[ci]] = s.upper - f.constant
end


function loadObj(cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.ScalarAffineFunction)
    for term in f.terms
        cbcModelFormat.obj[mapping.variable_mapping[term.variable_index]] = term.coefficient
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
    userOptimizer::MOI.AbstractOptimizer; copynames=false,
    warnattributes=false)


    nbCols = MOI.get(userOptimizer, MOI.NumberOfVariables())
    varIndex = MOI.get(userOptimizer, MOI.ListOfVariableIndices())
    for i in 1:nbCols
        cbcOptimizer.mapping.variable_mapping[varIndex[i]] = i
    end

    listOfConstraints = MOI.get(userOptimizer, MOI.ListOfConstraints())
    nbRows = 0;
    for (F,S) in listOfConstraints
        (F,S) in SUPPORTED_CONSTRAINTS || return MOI.CopyUnsupportedConstraint
        F == MOI.SingleVariable && continue
        ## Update constraint_mapping for (F,S) constraints
        ci = MOI.get(userOptimizer, MOI.ListOfConstraintIndices{F,S}())
        for i in 1:length(ci)
            cbcOptimizer.mapping.constraint_mapping[ci[i]] = nbRows + i
        end
        nbRows += MOI.get(userOptimizer, MOI.NumberOfConstraints{cons}())
    end

    cbcModelFormat = CbcModelFormat(nbRows, nbCols)

    for (F,S) in listOfConstraints
        for ci in MOI.get(userOptimizer, MOI.ListOfConstraintIndices{F,S}())
            f = MOI.get(userOptimizer, MOI.ConstraintFunction(), ci)
            s = MOI.get(userOptimizer,  MOI.ConstraintSet(), ci)
            S == MOI.Integer || loadConstraint(ci, cbcModelFormat, cbcOptimizer.mapping, f, s)
            if ((S == MOI.ZeroOne || S == MOI.Integer) && F == MOI.SingleVariable)
                setInteger(cbcOptimizer.inner, cbcOptimizer.mapping.variable_mapping[F.variable])
            end
        end
    end


    ## Copy objective function
    f = MOI.get(userOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    loadObj(cbcModelFormat, cbcOptimizer.mapping, f)


    loadProblem(cbcOptimizer.inner, cbcModelFormat.constraint_matrix, cbcModelFormat.col_lb, cbcModelFormat.col_ub, cbcModelFormat.obj, cbcModelFormat.row_lb, cbcModelFormat.row_ub)

    return MOI.CopySuccess
end


function MOI.optimize!(userOptimizer::CbcOptimizer)
    # Call solve function
    solve(userOptimizer.inner)
    return CbcCInterface.status(userOptimizer.inner)
end


# Set functions

function MOI.set!(userOptimizer::CbcOptimizer, object::MOI.ObjectiveSense, sense)
    sense == MOI.MinSense && setObjSense(userOptimizer.inner, 1)
    setObjSense(userOptimizer.inner, -1)
end


# Results-related getter functions

function MOI.get(userOptimizer::CbcOptimizer, object::MOI.ObjectiveValue)
    return getObjValue(userOptimizer.inner)
end

function MOI.get(userOptimizer::CbcOptimizer, object::MOI.ObjectiveBound)
    return getBestPossibleObjValue(userOptimizer.inner)
end



#
