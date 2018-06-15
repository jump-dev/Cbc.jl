export CbcOptimizer, copy!, optimize!, get


using MathOptInterface
using Cbc.CbcCInterface

const MOI = MathOptInterface
const CbcCI = CbcCInterface

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

const BOUND_CONSTRAINTS = [
    (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}),
    (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}),
    (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}),
    (MOI.SingleVariable, MOI.EqualTo{Float64}),
    (MOI.SingleVariable, MOI.LessThan{Float64}),
    (MOI.SingleVariable, MOI.GreaterThan{Float64}),
    (MOI.SingleVariable, MOI.Interval{Float64})
]

const INTEGRALITY_CONSTRAINTS = [
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
    constraint_matrix::Array{Float64,2}
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
    f::MOI.ScalarAffineFunction, s::MOI.EqualTo)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.constraint_mapping[ci],mapping.variable_mapping[term.variable_index]] = term.coefficient
    end
    cbcModelFormat.row_lb[mapping.constraint_mapping[ci]] = s.value - f.constant
    cbcModelFormat.row_ub[mapping.constraint_mapping[ci]] = s.value - f.constant
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.ScalarAffineFunction, s::MOI.GreaterThan)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.constraint_mapping[ci],mapping.variable_mapping[term.variable_index]] = term.coefficient
    end
    cbcModelFormat.row_lb[mapping.constraint_mapping[ci]] = s.lower - f.constant
end

function loadConstraint(ci::MOI.ConstraintIndex, cbcModelFormat::CbcModelFormat, mapping::VarAndConstrMapping,
    f::MOI.ScalarAffineFunction, s::MOI.LessThan)
    for term in f.terms
        cbcModelFormat.constraint_matrix[mapping.constraint_mapping[ci], mapping.variable_mapping[term.variable_index]] = term.coefficient
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

    zeroOneIndices = Vector{Int}(0)
    integerIndices = Vector{Int}(0)
    listOfConstraints = MOI.get(userOptimizer, MOI.ListOfConstraints())
    nbRows = 0;
    for (F,S) in listOfConstraints
        (F,S) in SUPPORTED_CONSTRAINTS || return MOI.CopyUnsupportedConstraint

        ci = MOI.get(userOptimizer, MOI.ListOfConstraintIndices{F,S}())

        if (F,S) in INTEGRALITY_CONSTRAINTS
            for i in 1:length(ci)
                f = MOI.get(userOptimizer, MOI.ConstraintFunction(), ci[i])
                S == MOI.ZeroOne && push!(zeroOneIndices, cbcOptimizer.mapping.variable_mapping[f.variable])
                S == MOI.Integer && push!(integerIndices, cbcOptimizer.mapping.variable_mapping[f.variable])
            end
        end

        ## Update constraint_mapping for (F,S) constraints
        F == MOI.SingleVariable && continue
        for i in 1:length(ci)
            cbcOptimizer.mapping.constraint_mapping[ci[i]] = nbRows + i
        end
        nbRows += MOI.get(userOptimizer, MOI.NumberOfConstraints{F,S}())

    end

    cbcModelFormat = CbcModelFormat(nbRows, nbCols)

    for (F,S) in listOfConstraints
        (F,S) in INTEGRALITY_CONSTRAINTS && continue
        for ci in MOI.get(userOptimizer, MOI.ListOfConstraintIndices{F,S}())
            f = MOI.get(userOptimizer, MOI.ConstraintFunction(), ci)
            s = MOI.get(userOptimizer,  MOI.ConstraintSet(), ci)
            loadConstraint(ci, cbcModelFormat, cbcOptimizer.mapping, f, s)
        end
    end

    ## Add bounds to binary columns if it makes sense
    for colIdx in zeroOneIndices
        if cbcModelFormat.col_lb[colIdx] < 0.0
            cbcModelFormat.col_lb[colIdx] = 0.0
        end
        if cbcModelFormat.col_ub[colIdx] > 1.0
            cbcModelFormat.col_ub[colIdx] = 1.0
        end
    end


    ## Copy objective function
    f = MOI.get(userOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    loadObj(cbcModelFormat, cbcOptimizer.mapping, f)

    # @show cbcModelFormat.constraint_matrix
    # @show cbcModelFormat.col_lb
    # @show cbcModelFormat.col_ub
    # @show cbcModelFormat.row_lb
    # @show cbcModelFormat.row_ub
    # @show cbcModelFormat.obj

    ## Load the problem to Cbc
    CbcCI.loadProblem(cbcOptimizer.inner, cbcModelFormat.constraint_matrix, cbcModelFormat.col_lb, cbcModelFormat.col_ub, cbcModelFormat.obj, cbcModelFormat.row_lb, cbcModelFormat.row_ub)

    ## Set integer variables
    for colIdx in vcat(integerIndices, zeroOneIndices)
        CbcCI.setInteger(cbcOptimizer.inner, colIdx-1)
    end

    return MOI.CopySuccess
end


function MOI.optimize!(cbcOptimizer::CbcOptimizer)
    # Call solve function
    CbcCI.solve(cbcOptimizer.inner)
    return CbcCI.status(cbcOptimizer.inner)
end


# Set functions

function MOI.set!(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveSense, sense)
    sense == MOI.MinSense && CbcCI.setObjSense(cbcOptimizer.inner, 1)
    sense == MOI.MaxSense && CbcCI.setObjSense(cbcOptimizer.inner, -1)
end


# Get functions

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.VariablePrimal, ref::MOI.VariableIndex)
    variablePrimals = CbcCI.getColSolution(cbcOptimizer.inner)
    return variablePrimals[cbcOptimizer.mapping.variable_mapping[ref]]
end

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveSense)
    CbcCI.getObjSense(cbcOptimizer.inner) == 1 && return MOI.MinSense
    CbcCI.getObjSense(cbcOptimizer.inner) == -1 && return MOI.MaxSense
end

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveValue)
    return CbcCI.getObjValue(cbcOptimizer.inner)
end

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.ObjectiveBound)
    return CbcCI.getBestPossibleObjValue(cbcOptimizer.inner)
end

function MOI.get(cbcOptimizer::CbcOptimizer, object::MOI.NodeCount)
    return CbcCI.getNodeCount(cbcOptimizer.inner)
end



#
