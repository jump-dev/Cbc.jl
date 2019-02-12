using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

import Cbc.CbcCInterface
const CbcCI = CbcCInterface

using Compat.SparseArrays

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::CbcCI.CbcModel
    objective_constant::Float64
    Optimizer() = new(CbcCI.CbcModel(), 0.0)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "COIN Branch-and-Cut (Cbc)"

function MOI.empty!(model::Optimizer)
    model.inner = CbcCI.CbcModel()
    model.objective_constant = 0.0
    CbcCI.setParameter(model.inner, "logLevel", "0")
    return
end

function MOI.is_empty(model::Optimizer)
    return CbcCI.getNumCols(model.inner) == 0 &&
        CbcCI.getNumRows(model.inner) == 0
end

MOI.add_variable(::Optimizer) = throw(MOI.AddVariableNotAllowed())

mutable struct CbcModelFormat
    objective_constant::Float64
    num_rows::Int
    num_cols::Int
    row_idx::Vector{Int}
    col_idx::Vector{Int}
    values::Vector{Float64}
    col_lb::Vector{Float64}
    col_ub::Vector{Float64}
    obj::Vector{Float64}
    row_lb::Vector{Float64}
    row_ub::Vector{Float64}
    binary::Vector{Int}
    integer::Vector{Int}
    function CbcModelFormat(num_rows::Integer, num_cols::Integer)
        obj = fill(0.0, num_cols)
        row_idx = Int[]
        col_idx = Int[]
        values = Float64[]
        col_lb = fill(-Inf, num_cols)
        col_ub = fill(Inf, num_cols)
        row_lb = fill(-Inf, num_rows)
        row_ub = fill(Inf, num_rows)
        # An `InexactError` might occur if `num_rows` or `num_cols` is too
        # large, e.g., if `num_cols isa Int64` and is larger than 2^31 on a
        # 32-bit hardware
        new(0.0, num_rows, num_cols, row_idx, col_idx, values, col_lb, col_ub,
            obj, row_lb, row_ub, Int[], Int[])
    end
end

###
### SingleVariable-in-{EqualTo, LessThan, GreaterThan, Interval, ZeroOne, Integer}
###

column_value(map::MOIU.IndexMap, f::MOI.VariableIndex) = map.varmap[f].value
column_value(map::MOIU.IndexMap, f::MOI.SingleVariable) = column_value(map, f.variable)

function MOI.supports_constraint(
        ::Optimizer, ::Type{MOI.SingleVariable}, ::Type{<:Union{
            MOI.EqualTo{Float64}, MOI.LessThan{Float64},
            MOI.GreaterThan{Float64}, MOI.Interval{Float64},
            MOI.ZeroOne, MOI.Integer}})
    return true
end

function load_constraint(
        ::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.SingleVariable, s::MOI.EqualTo)
    column = column_value(mapping, f)
    model.col_lb[column] = s.value
    model.col_ub[column] = s.value
    return
end

function load_constraint(
        ::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.SingleVariable, s::MOI.LessThan)
    column = column_value(mapping, f)
    model.col_ub[column] = s.upper
    return
end

function load_constraint(
        ::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.SingleVariable, s::MOI.GreaterThan)
    column = column_value(mapping, f)
    model.col_lb[column] = s.lower
    return
end

function load_constraint(
        ::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.SingleVariable, s::MOI.Interval)
    column = column_value(mapping, f)
    model.col_lb[column] = s.lower
    model.col_ub[column] = s.upper
    return
end

function load_constraint(
        ::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.SingleVariable, s::MOI.ZeroOne)
    push!(model.binary, column_value(mapping, f))
    return
end

function load_constraint(
        ::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.SingleVariable, s::MOI.Integer)
    push!(model.integer, column_value(mapping, f))
    return
end

###
### ScalarAffineFunction-in-{EqualTo, LessThan, GreaterThan, Interval}
###

function MOI.supports_constraint(
        ::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{<:Union{
            MOI.EqualTo{Float64}, MOI.LessThan{Float64},
            MOI.GreaterThan{Float64}, MOI.Interval{Float64}}})
    return true
end

function add_terms(
        model::CbcModelFormat, mapping::MOIU.IndexMap,
        ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S},
        func::MOI.ScalarAffineFunction{Float64}) where {S}
    for term in func.terms
        push!(model.row_idx, mapping.conmap[ci].value)
        push!(model.col_idx, column_value(mapping, term.variable_index))
        push!(model.values, term.coefficient)
    end
    return
end

function load_constraint(
        ci::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.ScalarAffineFunction, s::MOI.EqualTo)
    add_terms(model, mapping, ci, f)
    row = mapping.conmap[ci].value
    model.row_lb[row] = s.value - f.constant
    model.row_ub[row] = s.value - f.constant
    return
end

function load_constraint(
        ci::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.ScalarAffineFunction, s::MOI.GreaterThan)
    add_terms(model, mapping, ci, f)
    row = mapping.conmap[ci].value
    model.row_lb[row] = s.lower - f.constant
    return
end

function load_constraint(
        ci::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.ScalarAffineFunction, s::MOI.LessThan)
    add_terms(model, mapping, ci, f)
    row = mapping.conmap[ci].value
    model.row_ub[row] = s.upper - f.constant
    return
end

function load_constraint(
        ci::MOI.ConstraintIndex, model::CbcModelFormat,
        mapping::MOIU.IndexMap, f::MOI.ScalarAffineFunction, s::MOI.Interval)
    add_terms(model, mapping, ci, f)
    row = mapping.conmap[ci].value
    model.row_ub[row] = s.upper - f.constant
    model.row_lb[row] = s.lower - f.constant
    return
end


###
### ObjectiveSense
###

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense,
                 sense::MOI.OptimizationSense)
    if sense == MOI.MAX_SENSE
        CbcCI.setObjSense(model.inner, -1)
    else ## Other senses are set as minimization (cbc default)
        CbcCI.setObjSense(model.inner, 1)
    end
end

###
### ObjectiveFunction{ScalarAffineFunction}
###

function MOI.supports(
        ::Optimizer, ::MOI.ObjectiveFunction{
            <:Union{MOI.ScalarAffineFunction{Float64}, MOI.SingleVariable}})
    return true
end

function load_objective(
        model::CbcModelFormat, mapping::MOIU.IndexMap, f::MOI.ScalarAffineFunction)
    # We need to increment values of objective function with += to handle
    # cases like $x_1 + x_2 + x_1$. This is safe because objective function
    # is initialized with zeros in the constructor and `load_objective` only
    # gets called once.
    for term in f.terms
        column = column_value(mapping, term.variable_index)
        model.obj[column] += term.coefficient
    end
    model.objective_constant = f.constant
    return
end

function load_objective(
        model::CbcModelFormat, mapping::MOIU.IndexMap, f::MOI.SingleVariable)
    column = column_value(mapping, f)
    model.obj[column] = 1.0
    model.objective_constant = 0.0
    return
end

###
### This main copy_to function.
###

"""
    create_constraint_indices(src::MOI.ModelLike, mapping::MOIU.IndexMap)

Create a new set of constraint indices. Importantly:

 - The `.value` field of each `ConstraintIndex{ScalarAffineFunction, S}` is its
   row in the constraint matrix.
 - The `.value` field of each `ConstraintIndex{SingleVariable, S}` is the
   column of the associated variable in the constraint matrix.
"""
function create_constraint_indices(src::MOI.ModelLike, mapping::MOIU.IndexMap)
    num_rows = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        for index in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            if F == MOI.SingleVariable
                c_func = MOI.get(src, MOI.ConstraintFunction(), index)
                column = mapping.varmap[c_func.variable].value
                mapping.conmap[index] = MOI.ConstraintIndex{F, S}(column)
            else
                num_rows += 1
                mapping.conmap[index] = MOI.ConstraintIndex{F, S}(num_rows)
            end
        end
    end
    return num_rows
end

"""
    create_variable_indices(src::MOI.ModelLike, mapping::MOIU.IndexMap)

Create a new set of variable indices. Importantly, the `.value` field of each
`VariableIndex` is its column in the constraint matrix.
"""
function create_variable_indices(src::MOI.ModelLike, mapping::MOIU.IndexMap)
    var_index = MOI.get(src, MOI.ListOfVariableIndices())
    for (column, variable) in enumerate(var_index)
        mapping.varmap[variable] = MOI.VariableIndex(column)
    end
    return length(var_index)
end

function MOI.copy_to(cbc_dest::Optimizer, src::MOI.ModelLike;
                     copy_names=false)
    # Begin by checking that we support all the constraints. This means that
    # we don't have to do this later.
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if !(MOI.supports_constraint(cbc_dest, F, S))
            throw(MOI.UnsupportedConstraint{F, S}(
                "Cbc.Optimizer does not support constraints of type $F -in- $S."
            ))
        end
    end

    # A mapping between the indices of `src` and the ones we create.
    mapping = MOIU.IndexMap()

    num_cols = create_variable_indices(src, mapping)
    num_rows = create_constraint_indices(src, mapping)

    # Create a new temporary storage instance.
    tmp_model = CbcModelFormat(num_rows, num_cols)

    # Copy the constraints out of `src` into `tmp_model`.
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        for index in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            load_constraint(
                index,
                tmp_model,
                mapping,
                MOI.get(src, MOI.ConstraintFunction(), index),
                MOI.get(src,  MOI.ConstraintSet(), index)
            )
        end
    end

    # Since Cbc doesn't have an explicit binary variable, we need to add [0, 1]
    # bounds and make it integer (which is done at the end of this function).
    for column in tmp_model.binary
        tmp_model.col_lb[column] = max(tmp_model.col_lb[column], 0.0)
        tmp_model.col_ub[column] = min(tmp_model.col_ub[column], 1.0)
    end

    # Copy the objective function.
    objective_function_type = MOI.get(src, MOI.ObjectiveFunctionType())
    if !MOI.supports(cbc_dest, MOI.ObjectiveFunction{objective_function_type}())
        error("Objective function type $(objective_function_type) not supported.")
    end
    objective_function = MOI.get(
        src, MOI.ObjectiveFunction{objective_function_type}())
    load_objective(tmp_model, mapping, objective_function)

    sense = MOI.get(src, MOI.ObjectiveSense())
    MOI.set(cbc_dest, MOI.ObjectiveSense(), sense)

    # Load the problem into Cbc.
    CbcCI.loadProblem(
        cbc_dest.inner,
        sparse(tmp_model.row_idx, tmp_model.col_idx, tmp_model.values,
               tmp_model.num_rows, tmp_model.num_cols),
        tmp_model.col_lb, tmp_model.col_ub,
        tmp_model.obj,
        tmp_model.row_lb, tmp_model.row_ub
    )

    cbc_dest.objective_constant = tmp_model.objective_constant

    # Set the integer variables.
    for idx in tmp_model.integer
        CbcCI.setInteger(cbc_dest.inner, idx - 1)
    end

    # Set the binary variables.
    for idx in tmp_model.binary
        CbcCI.setInteger(cbc_dest.inner, idx - 1)
    end

    return MOIU.IndexMap(mapping.varmap, mapping.conmap)
end

###
### Optimize and post-optimize functions
###

function MOI.optimize!(model::Optimizer)
    CbcCI.solve(model.inner)
    return
end

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return CbcCI.getNumCols(model.inner)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    return CbcCI.getBestPossibleObjValue(model.inner) + model.objective_constant
end

function MOI.get(model::Optimizer, ::MOI.NodeCount)
    return CbcCI.getNodeCount(model.inner)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
    return CbcCI.getObjValue(model.inner) + model.objective_constant
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimal,
                 ref::MOI.VariableIndex)
    primal_solution = CbcCI.unsafe_getColSolution(model.inner)
    return primal_solution[ref.value]
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimal,
                 ref::Vector{MOI.VariableIndex})
    primal_solution = CbcCI.unsafe_getColSolution(model.inner)
    return [primal_solution[vi.value] for vi in ref]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any})
    primal_solution = CbcCI.unsafe_getRowActivity(model.inner)
    return primal_solution[ci.value]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, <:Any})
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(ci.value))
end

function MOI.get(model::Optimizer, object::MOI.ResultCount)
    if CbcCI.isProvenInfeasible(model.inner) ||
            CbcCI.isContinuousUnbounded(model.inner) ||
            CbcCI.isAbandoned(model.inner) ||
            CbcCI.getObjValue(model.inner) >= 1e300
        return 0
    end
    return 1
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    sense = CbcCI.getObjSense(model.inner)
    if sense == 1
        return MOI.MIN_SENSE
    elseif sense == -1
        return MOI.MAX_SENSE
    else
        error("Objective sense $(sense) not recognized.")
    end
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if CbcCI.isProvenInfeasible(model.inner)
        return MOI.INFEASIBLE
    elseif CbcCI.isContinuousUnbounded(model.inner)
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif CbcCI.isNodeLimitReached(model.inner)
        return MOI.NODE_LIMIT
    elseif CbcCI.isSecondsLimitReached(model.inner)
        return MOI.TIME_LIMIT
    elseif CbcCI.isSolutionLimitReached(model.inner)
        return MOI.SOLUTION_LIMIT
    elseif CbcCI.isProvenOptimal(model.inner) ||
            CbcCI.isInitialSolveProvenOptimal(model.inner) ||
            MOI.get(model, MOI.ResultCount()) == 1
        return MOI.OPTIMAL
    elseif CbcCI.isAbandoned(model.inner)
        return MOI.INTERRUPTED
    elseif CbcCI.optimizeNotCalled(model.inner)
        return MOI.OPTIMIZE_NOT_CALLED
    else
        error("Internal error: Unrecognized solution status:",
              " status = $(status(model.inner)),",
              " secondaryStatus = $(secondaryStatus(model.inner))")
    end
end

function MOI.get(model::Optimizer, ::MOI.PrimalStatus)
    if CbcCI.isProvenOptimal(model.inner) ||
            CbcCI.isInitialSolveProvenOptimal(model.inner)
        return MOI.FEASIBLE_POINT
    elseif CbcCI.isProvenInfeasible(model.inner)
        return MOI.INFEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(::Optimizer, ::MOI.DualStatus)
    return MOI.NO_SOLUTION
end
