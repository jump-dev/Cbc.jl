import MathOptInterface
import SparseArrays

const MOI = MathOptInterface

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Ptr{Cvoid}
    silent::Bool
    params::Dict{String, String}
    variable_start::Dict{MOI.VariableIndex, Float64}
    objective_constant::Float64
    solve_time::Float64

    """
        Optimizer(; kwargs...)

    Create a new Cbc Optimizer.
    """
    function Optimizer(; kwargs...)
        model = new(
            Cbc_newModel(),
            false,
            Dict{String, String}(),
            Dict{MOI.VariableIndex, Float64}(),
            0.0,
            0.0,
        )
        for (key, value) in kwargs
            MOI.set(model, MOI.RawParameter(key), value)
        end
        finalizer(model) do m
            Cbc_deleteModel(m.inner)
        end
        return model
    end
end

function MOI.supports(::Optimizer, ::MOI.RawParameter)
    # TODO(odow): There is no programatical way throught the C API to check if a
    # parameter name (or value) is valid. Fix this upstream.
    return true
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    return MOI.set(model, param, string(value))
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value::String)
    if !MOI.supports(model, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    name = string(param.name)
    model.params[name] = value
    if !(model.silent && name == "logLevel")
        Cbc_setParameter(model.inner, name, value)
    end
    return
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    # TODO: This gives a poor error message if the name of the parameter is invalid.
    return model.params[string(param.name)]
end

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    if value == model.silent
        return
    end
    log_level = value ? "0" : get(model.params, "logLevel", "1")
    Cbc_setParameter(model.inner, "logLevel", log_level)
    model.silent = value
    return
end
MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value)
    if value === nothing
        delete!(model.params, "seconds")
        Cbc_setParameter(model.inner, "seconds", "??")
    else
        MOI.set(model, MOI.RawParameter("seconds"), value)
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    value = get(model.params, "seconds", nothing)
    return value === nothing ? value : parse(Float64, value)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "COIN Branch-and-Cut (Cbc)"

function MOI.empty!(model::Optimizer)
    Cbc_deleteModel(model.inner)
    model.inner = Cbc_newModel()
    model.objective_constant = 0.0
    model.solve_time = 0.0
    for (name, value) in model.params
        Cbc_setParameter(model.inner, name, value)
    end
    if model.silent
        Cbc_setParameter(model.inner, "logLevel", "0")
    end
    empty!(model.variable_start)
    return
end

function MOI.is_empty(model::Optimizer)
    return Cbc_getNumCols(model.inner) == 0 && Cbc_getNumRows(model.inner) == 0
end

mutable struct CbcModelFormat
    num_rows::Cint
    num_cols::Cint
    # (row_idx, col_idx, values) are the sparse elements of the constraint matrix.
    row_idx::Vector{Cint}
    col_idx::Vector{Cint}
    values::Vector{Float64}
    # Constraint bounds.
    row_lb::Vector{Float64}
    row_ub::Vector{Float64}
    # Variable bounds.
    col_lb::Vector{Float64}
    col_ub::Vector{Float64}
    # Columns that are binary or integer
    binary::Vector{Int}
    integer::Vector{Int}
    # SOS1 constraints
    sos1_starts::Vector{Cint}
    sos1_indices::Vector{Cint}
    sos1_weights::Vector{Float64}
    # SOS2 constraints
    sos2_starts::Vector{Cint}
    sos2_indices::Vector{Cint}
    sos2_weights::Vector{Float64}
    # Objective coefficients.
    obj::Vector{Float64}
    objective_constant::Float64

    # An `InexactError` might occur if `num_rows` or `num_cols` is too
    # large, e.g., if `num_cols isa Int64` and is larger than 2^31 on
    # 32-bit hardware.
    function CbcModelFormat(num_rows::Cint, num_cols::Cint)
        return new(
            num_rows,
            num_cols,
            # Constraint matrix.
            Cint[],
            Cint[],
            Float64[],
            # Row lower/upper bounds.
            fill(-Inf, num_rows),
            fill(Inf, num_rows),
            # Column lower/upper bounds.
            fill(-Inf, num_cols),
            fill(Inf, num_cols),
            # Binary/Integer columns.
            Cint[],
            Cint[],
            # SOSI constraints.
            Cint[],
            Cint[],
            Float64[],
            # SOSII constraints.
            Cint[],
            Cint[],
            Float64[],
            # Objective vector and offset.
            fill(0.0, num_cols),
            0.0,
        )
    end
end

###
### SingleVariable-in-{EqualTo,LessThan,GreaterThan,Interval,ZeroOne,Integer}
###

function column_value(map::MOI.Utilities.IndexMap, index::MOI.VariableIndex)
    return map.varmap[index].value
end

function column_value(map::MOI.Utilities.IndexMap, func::MOI.SingleVariable)
    return column_value(map, func.variable)
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.SingleVariable},
    ::Type{
         <:Union{
            MOI.EqualTo{Float64},
            MOI.LessThan{Float64},
            MOI.GreaterThan{Float64},
            MOI.Interval{Float64},
            MOI.ZeroOne,
            MOI.Integer,
        }
    }
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}}
)
    return true
end

function load_constraint(
    ::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.SingleVariable,
    set::MOI.EqualTo
)
    column = column_value(mapping, func)
    model.col_lb[column] = set.value
    model.col_ub[column] = set.value
    return
end

function load_constraint(
    ::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.SingleVariable,
    set::MOI.LessThan
)
    column = column_value(mapping, func)
    model.col_ub[column] = set.upper
    return
end

function load_constraint(
    ::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.SingleVariable,
    set::MOI.GreaterThan
)
    column = column_value(mapping, func)
    model.col_lb[column] = set.lower
    return
end

function load_constraint(
    ::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.SingleVariable,
    set::MOI.Interval
)
    column = column_value(mapping, func)
    model.col_lb[column] = set.lower
    model.col_ub[column] = set.upper
    return
end

function load_constraint(
    ::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.SingleVariable,
    ::MOI.ZeroOne
)
    push!(model.binary, column_value(mapping, func) - 1)
    return
end

function load_constraint(
    ::MOI.ConstraintIndex,
     model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.SingleVariable,
    ::MOI.Integer
)
    push!(model.integer, column_value(mapping, func) - 1)
    return
end

function load_constraints(
    model::CbcModelFormat,
    src::MOI.ModelLike,
    mapping::MOI.Utilities.IndexMap,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet}
)
    for index in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        load_constraint(
            index,
            model,
            mapping,
            MOI.get(src, MOI.ConstraintFunction(), index),
            MOI.get(src,  MOI.ConstraintSet(), index)
        )
    end
    return
end

###
### ScalarAffineFunction-in-{EqualTo, LessThan, GreaterThan, Interval}
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{
        <:Union{
            MOI.EqualTo{Float64},
            MOI.LessThan{Float64},
            MOI.GreaterThan{Float64},
            MOI.Interval{Float64},
        }
    }
)
    return true
end

function add_terms(
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    index::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S},
    func::MOI.ScalarAffineFunction{Float64},
) where {S}
    for term in func.terms
        push!(model.row_idx, mapping.conmap[index].value)
        push!(model.col_idx, column_value(mapping, term.variable_index))
        push!(model.values, term.coefficient)
    end
    return
end

function load_constraint(
    index::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.ScalarAffineFunction,
    set::MOI.EqualTo,
)
    add_terms(model, mapping, index, func)
    row = mapping.conmap[index].value
    model.row_lb[row] = set.value - func.constant
    model.row_ub[row] = set.value - func.constant
    return
end

function load_constraint(
    index::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.ScalarAffineFunction,
    set::MOI.GreaterThan,
)
    add_terms(model, mapping, index, func)
    row = mapping.conmap[index].value
    model.row_lb[row] = set.lower - func.constant
    return
end

function load_constraint(
    index::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.ScalarAffineFunction,
    set::MOI.LessThan,
)
    add_terms(model, mapping, index, func)
    row = mapping.conmap[index].value
    model.row_ub[row] = set.upper - func.constant
    return
end

function load_constraint(
    index::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.ScalarAffineFunction,
    set::MOI.Interval,
)
    add_terms(model, mapping, index, func)
    row = mapping.conmap[index].value
    model.row_ub[row] = set.upper - func.constant
    model.row_lb[row] = set.lower - func.constant
    return
end

function load_constraint(
    index::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.VectorOfVariables,
    set::MOI.SOS1{Float64},
)
    push!(model.sos1_starts, Cint(length(model.sos1_weights)))
    append!(model.sos1_weights, set.weights)
    for v in func.variables
        push!(model.sos1_indices, Cint(v.value - 1))
    end
    return
end

function load_constraint(
    index::MOI.ConstraintIndex,
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    func::MOI.VectorOfVariables,
    set::MOI.SOS2{Float64},
)
    push!(model.sos2_starts, Cint(length(model.sos2_weights)))
    append!(model.sos2_weights, set.weights)
    for v in func.variables
        push!(model.sos2_indices, Cint(v.value - 1))
    end
    return
end

###
### ObjectiveSense
###

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense
)
    if sense == MOI.MAX_SENSE
        Cbc_setObjSense(model.inner, -1.0)
    elseif sense == MOI.MIN_SENSE
        Cbc_setObjSense(model.inner, 1.0)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        Cbc_setObjSense(model.inner, 0.0)
    end
    return
end

###
### ObjectiveFunction{ScalarAffineFunction}
###

function MOI.supports(
    ::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
)
    return true
end

function load_objective(
    model::CbcModelFormat,
    mapping::MOI.Utilities.IndexMap,
    dest::Optimizer,
    src::MOI.ModelLike,
)
    F = MOI.get(src, MOI.ObjectiveFunctionType())
    if !MOI.supports(dest, MOI.ObjectiveFunction{F}())
        error("Objective function type $(F) not supported.")
    end
    f = MOI.get(src, MOI.ObjectiveFunction{F}())
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

###
### Variable starting values
###

function MOI.supports(
    ::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}
)
    return true
end

function MOI.set(
    model::Optimizer, ::MOI.VariablePrimalStart, x::MOI.VariableIndex, ::Nothing
)
    delete!(model.variable_start, x)
    return
end

function MOI.set(
    model::Optimizer, ::MOI.VariablePrimalStart, x::MOI.VariableIndex, value
)
    model.variable_start[x] = value
    return
end

function MOI.get(
    model::Optimizer, ::MOI.VariablePrimalStart, x::MOI.VariableIndex
)
    return get(model.variable_start, x, nothing)
end

###
### This main copy_to function.
###

function create_constraint_indices_for_types(
    src::MOI.ModelLike,
    mapping::MOI.Utilities.IndexMap,
    ::Type{MOI.SingleVariable},
    S::Type{<:MOI.AbstractSet},
    num_rows::Int,
)
    for index in MOI.get(
        src, MOI.ListOfConstraintIndices{MOI.SingleVariable, S}()
    )
        f = MOI.get(src, MOI.ConstraintFunction(), index)
        i = mapping.varmap[f.variable].value
        mapping.conmap[index] = MOI.ConstraintIndex{MOI.SingleVariable, S}(i)
    end
    return num_rows
end

function create_constraint_indices_for_types(
    src::MOI.ModelLike,
    mapping::MOI.Utilities.IndexMap,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
    num_rows::Int,
)
    for index in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        num_rows += 1
        mapping.conmap[index] = MOI.ConstraintIndex{F, S}(num_rows)
    end
    return num_rows
end

function create_constraint_indices(
    dest::Optimizer, src::MOI.ModelLike, mapping::MOI.Utilities.IndexMap
)
    n = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if !(MOI.supports_constraint(dest, F, S))
            throw(MOI.UnsupportedConstraint{F, S}(
                "Cbc.Optimizer does not support constraints of type $F -in- $S."
            ))
        end
        # The type of `F` and `S` is not type-stable, so we use a function
        # barrier (`create_constraint_indices_for_types`) to improve performance.
        n = create_constraint_indices_for_types(src, mapping, F, S, n)
    end
    return Cint(n)
end

function create_variable_indices(
    src::MOI.ModelLike, mapping::MOI.Utilities.IndexMap
)
    for (i, x) in enumerate(MOI.get(src, MOI.ListOfVariableIndices()))
        mapping.varmap[x] = MOI.VariableIndex(i)
    end
    return Cint(length(mapping.varmap))
end

function MOI.copy_to(
    cbc_dest::Optimizer, src::MOI.ModelLike; copy_names::Bool = false
)
    @assert MOI.is_empty(cbc_dest)
    mapping = MOI.Utilities.IndexMap()
    num_cols = create_variable_indices(src, mapping)
    num_rows = create_constraint_indices(cbc_dest, src, mapping)
    tmp_model = CbcModelFormat(num_rows, num_cols)

    load_objective(tmp_model, mapping, cbc_dest, src)

    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        # The type of `F` and `S` is not type-stable, so we use a function
        # barrier (`load_constraints`) to improve performance.
        load_constraints(tmp_model, src, mapping, F, S)
    end

    # Since Cbc doesn't have an explicit binary variable, we need to add [0, 1]
    # bounds and make it integer (which is done at the end of this function).
    for column in tmp_model.binary
        tmp_model.col_lb[column + 1] = max(tmp_model.col_lb[column + 1], 0.0)
        tmp_model.col_ub[column + 1] = min(tmp_model.col_ub[column + 1], 1.0)
    end

    A = SparseArrays.sparse(
        tmp_model.row_idx,
        tmp_model.col_idx,
        tmp_model.values,
        tmp_model.num_rows,
        tmp_model.num_cols,
    )

    Cbc_loadProblem(
        cbc_dest.inner,
        tmp_model.num_cols,
        tmp_model.num_rows,
        A.colptr .- Cint(1),
        A.rowval .- Cint(1),
        A.nzval,
        tmp_model.col_lb,
        tmp_model.col_ub,
        tmp_model.obj,
        tmp_model.row_lb,
        tmp_model.row_ub,
    )

    MOI.Utilities.pass_attributes(
        cbc_dest,
        src,
        copy_names,
        mapping,
        MOI.get(src, MOI.ListOfVariableIndices()),
    )

    for attr in MOI.get(src, MOI.ListOfModelAttributesSet())
        if attr isa MOI.ObjectiveFunction
            continue # Already copied
        elseif !copy_names && attr isa MOI.Name
            continue
        end
        value = MOI.get(src, attr)
        if value !== nothing
            mapped_value = MOI.Utilities.map_indices(mapping, value)
            MOI.set(cbc_dest, attr, mapped_value)
        end
    end

    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        MOI.Utilities.pass_attributes(
            cbc_dest, src, copy_names, mapping, cis_src
        )
    end

    cbc_dest.objective_constant = tmp_model.objective_constant

    Cbc_setInteger.(cbc_dest.inner, tmp_model.integer)
    Cbc_setInteger.(cbc_dest.inner, tmp_model.binary)

    if length(tmp_model.sos1_starts) > 0
        push!(tmp_model.sos1_starts, Cint(length(tmp_model.sos1_weights)))
        Cbc_addSOS(
            cbc_dest.inner,
            length(tmp_model.sos1_starts) - 1,
            tmp_model.sos1_starts,
            tmp_model.sos1_indices,
            tmp_model.sos1_weights,
            Cint(1),
        )
    end
    if length(tmp_model.sos2_starts) > 0
        push!(tmp_model.sos2_starts, Cint(length(tmp_model.sos2_weights)))
        Cbc_addSOS(
            cbc_dest.inner,
            length(tmp_model.sos2_starts) - 1,
            tmp_model.sos2_starts,
            tmp_model.sos2_indices,
            tmp_model.sos2_weights,
            Cint(2),
        )
    end
    return mapping
end

###
### Optimize and post-optimize functions
###

function _unsafe_wrap_cbc_array(
    model::Ptr{Cvoid},
    f::Function,
    n::Integer,
    indices;
    own::Bool = false
)
    p = f(model)
    if p == C_NULL
        return map(x -> NaN, indices)
    end
    x = unsafe_wrap(Array, p, (n,); own = own)
    return x[indices]
end

function MOI.optimize!(model::Optimizer)
    if !isempty(model.variable_start)
        columns = fill(Cint(0), length(model.variable_start))
        values = fill(0.0, length(model.variable_start))
        for (i, (variable, value)) in enumerate(model.variable_start)
            columns[i] = Cint(variable.value - 1)
            values[i] = value
        end
        Cbc_setMIPStartI(model.inner, length(columns), columns, values)
    end
    t = time()
    Cbc_solve(model.inner)
    model.solve_time = time() - t
    return
end

function MOI.get(model::Optimizer, ::MOI.SolveTime)
    return model.solve_time
end

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return Cbc_getNumCols(model.inner)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    return Cbc_getBestPossibleObjValue(model.inner) + model.objective_constant
end

function MOI.get(model::Optimizer, ::MOI.NodeCount)
    return Cbc_getNodeCount(model.inner)
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return Cbc_getObjValue(model.inner) + model.objective_constant
end

function MOI.get(
    model::Optimizer, attr::MOI.VariablePrimal, x::MOI.VariableIndex
)
    MOI.check_result_index_bounds(model, attr)
    return _unsafe_wrap_cbc_array(
        model.inner,
        Cbc_getColSolution,
        Cbc_getNumCols(model.inner),
        x.value
    )
end

function MOI.get(
    model::Optimizer, attr::MOI.VariablePrimal, x::Vector{MOI.VariableIndex}
)
    MOI.check_result_index_bounds(model, attr)
    return _unsafe_wrap_cbc_array(
        model.inner,
        Cbc_getColSolution,
        Cbc_getNumCols(model.inner),
        [xi.value for xi in x],
    )
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    index::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    MOI.check_result_index_bounds(model, attr)
    return _unsafe_wrap_cbc_array(
        model.inner,
        Cbc_getRowActivity,
        Cbc_getNumRows(model.inner),
        index.value
    )
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    index::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}
)
    MOI.check_result_index_bounds(model, attr)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(index.value))
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    sense = Cbc_getObjSense(model.inner)
    if sense == 1.0
        return MOI.MIN_SENSE
    elseif sense == -1.0
        return MOI.MAX_SENSE
    else
        @assert sense == 0.0
        return MOI.FEASIBILITY_SENSE
    end
end

struct Status <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::Status) = true

struct SecondaryStatus <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::SecondaryStatus) = true

MOI.get(model::Optimizer, ::Status) = Cbc_status(model.inner)
MOI.get(model::Optimizer, ::SecondaryStatus) = Cbc_secondaryStatus(model.inner)

const _STATUS = Dict{Cint, String}(
    Cint(-1) => "before branchAndBound",
    Cint(0) => "finished - check isProvenOptimal or isProvenInfeasible to see if solution found (or check value of best solution)",
    Cint(1) => "stopped - on maxnodes, maxsols, maxtime",
    Cint(2) => "execution abandoned due to numerical dificulties",
    Cint(5) => "user programmed interruption",
)

const _SECONDARY_STATUS = Dict{Cint, String}(
    Cint(-1) => "unset (status_ will also be -1)",
    Cint(0) => "search completed with solution",
    Cint(1) => "linear relaxation not feasible (or worse than cutoff)",
    Cint(2) => "stopped on gap",
    Cint(3) => "stopped on nodes",
    Cint(4) => "stopped on time",
    Cint(5) => "stopped on user event",
    Cint(6) => "stopped on solutions",
    Cint(7) => "linear relaxation unbounded",
    Cint(8) => "stopped on iteration limit",
)

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    return """
    Cbc_status          = $(_STATUS[Cbc_status(model.inner)])
    Cbc_secondaryStatus = $(_SECONDARY_STATUS[Cbc_secondaryStatus(model.inner)])
    """
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    return Cbc_bestSolution(model.inner) != C_NULL ? 1 : 0
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    status = Cbc_status(model.inner)
    if status == -1
        return MOI.OPTIMIZE_NOT_CALLED
    elseif status == 0
        if Cbc_isProvenOptimal(model.inner) != 0
            return MOI.OPTIMAL
        elseif Cbc_isProvenInfeasible(model.inner) != 0
            return MOI.INFEASIBLE
        elseif Cbc_isContinuousUnbounded(model.inner) != 0
            return MOI.INFEASIBLE_OR_UNBOUNDED
        else
            return MOI.OTHER_ERROR
        end
    elseif status == 1
        if Cbc_isNodeLimitReached(model.inner) != 0
            return MOI.NODE_LIMIT
        elseif Cbc_isSecondsLimitReached(model.inner) != 0
            return MOI.TIME_LIMIT
        elseif Cbc_isSolutionLimitReached(model.inner) != 0
            return MOI.SOLUTION_LIMIT
        else
            return MOI.OTHER_LIMIT
        end
    elseif status == 2
        return MOI.NUMERICAL_ERROR
    else
        @assert status == 5
        return MOI.INTERRUPTED
    end
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    end
    if MOI.get(model, MOI.ResultCount()) > 0
        return MOI.FEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(::Optimizer, ::MOI.DualStatus)
    return MOI.NO_SOLUTION
end
