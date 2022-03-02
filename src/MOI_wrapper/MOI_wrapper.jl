const MOI = MathOptInterface

MOI.Utilities.@product_of_sets(
    _LPProductOfSets,
    MOI.EqualTo{T},
    MOI.LessThan{T},
    MOI.GreaterThan{T},
    MOI.Interval{T},
)

MOI.Utilities.@struct_of_constraints_by_set_types(
    _CbcConstraints,
    Union{MOI.EqualTo{T},MOI.LessThan{T},MOI.GreaterThan{T},MOI.Interval{T}},
    MOI.SOS1{T},
    MOI.SOS2{T},
)

const OptimizerCache = MOI.Utilities.GenericModel{
    Float64,
    MOI.Utilities.ObjectiveContainer{Float64},
    MOI.Utilities.VariablesContainer{Float64},
    _CbcConstraints{Float64}{
        MOI.Utilities.MatrixOfConstraints{
            Float64,
            MOI.Utilities.MutableSparseMatrixCSC{
                Float64,
                Cint,
                MOI.Utilities.ZeroBasedIndexing,
            },
            MOI.Utilities.Hyperrectangle{Float64},
            _LPProductOfSets{Float64},
        },
        MOI.Utilities.VectorOfConstraints{
            MOI.VectorOfVariables,
            MOI.SOS1{Float64},
        },
        MOI.Utilities.VectorOfConstraints{
            MOI.VectorOfVariables,
            MOI.SOS2{Float64},
        },
    },
}

"""
    Optimizer()

Create a new Cbc Optimizer.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Ptr{Cvoid}
    silent::Bool
    params::Dict{String,String}
    variable_start::Dict{MOI.VariableIndex,Float64}
    objective_constant::Float64
    solve_time::Float64
    termination_status::Cint
    has_solution::Bool
    variable_primal::Union{Nothing,Vector{Float64}}
    constraint_primal::Union{Nothing,Vector{Float64}}

    function Optimizer()
        model = new(
            Cbc_newModel(),
            false,
            Dict{String,String}(),
            Dict{MOI.VariableIndex,Float64}(),
            0.0,
            0.0,
            Cint(-1),
            false,
            nothing,
            nothing,
        )
        finalizer(Cbc_deleteModel, model)
        return model
    end
end

Base.cconvert(::Type{Ptr{Cvoid}}, model::Optimizer) = model
Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Optimizer) = model.inner

function MOI.default_cache(::Optimizer, ::Type{Float64})
    return MOI.Utilities.UniversalFallback(OptimizerCache())
end

function MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute)
    # TODO(odow): There is no programatical way throught the C API to check if a
    # parameter name (or value) is valid. Fix this upstream.
    return true
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return MOI.set(model, param, string(value))
end

function MOI.set(
    model::Optimizer,
    param::MOI.RawOptimizerAttribute,
    value::String,
)
    if !MOI.supports(model, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    model.params[param.name] = value
    if name == "threads" && Sys.iswindows()
        @warn(
            "Ignoring threads parameter due to known bugs in Cbc.jl. Read " *
            "https://github.com/jump-dev/Cbc.jl/issues/186 for more details.",
        )
        return
    end
    if !(model.silent && name == "logLevel")
        Cbc_setParameter(model, name, value)
    end
    return
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    # TODO: This gives a poor error message if the name of the parameter is
    # invalid.
    return model.params[param.name]
end

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    if value == model.silent
        return
    end
    log_level = value ? "0" : get(model.params, "logLevel", "1")
    Cbc_setParameter(model, "logLevel", log_level)
    model.silent = value
    return
end
MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value)
    if value === nothing
        delete!(model.params, "seconds")
        Cbc_setParameter(model, "seconds", "??")
    else
        MOI.set(model, MOI.RawOptimizerAttribute("seconds"), value)
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    value = get(model.params, "seconds", nothing)
    return value === nothing ? value : parse(Float64, value)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "COIN Branch-and-Cut (Cbc)"

MOI.get(::Optimizer, ::MOI.SolverVersion) = unsafe_string(Cbc_getVersion())

function MOI.empty!(model::Optimizer)
    Cbc_deleteModel(model)
    model.inner = Cbc_newModel()
    model.objective_constant = 0.0
    model.termination_status = Cint(-1)
    model.solve_time = 0.0
    for (name, value) in model.params
        MOI.set(model, MOI.RawOptimizerAttribute(name), value)
    end
    if model.silent
        Cbc_setParameter(model, "logLevel", "0")
    end
    empty!(model.variable_start)
    model.has_solution = false
    model.variable_primal = nothing
    model.constraint_primal = nothing
    return
end

function MOI.is_empty(model::Optimizer)
    return Cbc_getNumCols(model) == 0 && Cbc_getNumRows(model) == 0
end

function _index_map(
    src::OptimizerCache,
    ::MOI.IndexMap,
    ci::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction},
)
    return MOI.Utilities.rows(src.constraints.moi_equalto, ci)
end

function _index_map(
    ::OptimizerCache,
    index_map::MOI.IndexMap,
    ci::MOI.ConstraintIndex{<:MOI.VariableIndex},
)
    return index_map[MOI.VariableIndex(ci.value)].value
end

function _index_map(
    ::OptimizerCache,
    ::MOI.IndexMap,
    ci::MOI.ConstraintIndex{<:MOI.VectorOfVariables},
)
    return ci.value
end

function _index_map(
    src::OptimizerCache,
    index_map::MOI.IndexMap,
    ::Type{F},
    ::Type{S},
) where {F,S}
    inner = index_map.con_map[F, S]
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        inner[ci] = MOI.ConstraintIndex{F,S}(_index_map(src, index_map, ci))
    end
    return
end

"""
    _index_map(src::OptimizerCache)

Create an `IndexMap` mapping the variables and constraints in `OptimizerCache`
to their corresponding 1-based columns and rows.
"""
function _index_map(src::OptimizerCache)
    index_map = MOI.IndexMap()
    for (i, x) in enumerate(MOI.get(src, MOI.ListOfVariableIndices()))
        index_map[x] = MOI.VariableIndex(i)
    end
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        _index_map(src, index_map, F, S)
    end
    return index_map
end

function _constraint_matrix(constraints, n)
    @assert n == constraints.coefficients.n
    return (
        m = constraints.coefficients.m,
        n = constraints.coefficients.n,
        colptr = constraints.coefficients.colptr,
        rowval = constraints.coefficients.rowval,
        nzval = constraints.coefficients.nzval,
        lower = constraints.constants.lower,
        upper = constraints.constants.upper,
    )
end

function _constraint_matrix(::Nothing, n)
    return (
        m = Cint(0),
        n = Cint(n),
        colptr = fill(Cint(0), n + 1),
        rowval = Cint[],
        nzval = Float64[],
        lower = Float64[],
        upper = Float64[],
    )
end

function MOI.copy_to(dest::Optimizer, src::OptimizerCache)
    matrix = _constraint_matrix(
        src.constraints.moi_equalto,
        MOI.get(src, MOI.NumberOfVariables()),
    )
    c = zeros(matrix.n)
    obj =
        MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    for term in obj.terms
        c[term.variable.value] += term.coefficient
    end
    dest.objective_constant = obj.constant
    zeroone_attr = MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.ZeroOne}()
    binaries = Cint[Cint(ci.value - 1) for ci in MOI.get(src, zeroone_attr)]
    variable_lower = copy(src.variables.lower)
    variable_upper = copy(src.variables.upper)
    for b in binaries
        variable_lower[b+1] = max(variable_lower[b+1], 0.0)
        variable_upper[b+1] = min(variable_upper[b+1], 1.0)
    end
    Cbc_loadProblem(
        dest,
        matrix.n,
        matrix.m,
        matrix.colptr,
        matrix.rowval,
        matrix.nzval,
        variable_lower,
        variable_upper,
        c,
        matrix.lower,
        matrix.upper,
    )
    sense = MOI.get(src, MOI.ObjectiveSense())
    if sense == MOI.MIN_SENSE
        Cbc_setObjSense(dest, 1)
    elseif sense == MOI.MAX_SENSE
        Cbc_setObjSense(dest, -1)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        Cbc_setObjSense(dest, 0)
    end
    Cbc_setInteger.(dest, binaries)
    attr = MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.Integer}()
    for ci in MOI.get(src, attr)
        Cbc_setInteger(dest, Cint(ci.value - 1))
    end
    any_sos = false
    for (S, type) in ((MOI.SOS1{Float64}, 1), (MOI.SOS2{Float64}, 2))
        starts, indices, weights = Cint[], Cint[], Float64[]
        attr = MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S}()
        for ci in MOI.get(src, attr)
            any_sos = true
            push!(starts, Cint(length(weights)))
            f = MOI.get(src, MOI.ConstraintFunction(), ci)
            for x in f.variables
                push!(indices, Cint(x.value - 1))
            end
            s = MOI.get(src, MOI.ConstraintSet(), ci)
            append!(weights, s.weights)
        end
        N = Cint(length(starts))
        if N > 0
            push!(starts, length(weights))
            Cbc_addSOS(dest, N, starts, indices, weights, Cint(type))
        end
    end
    if any_sos && Cbc_getNumIntegers(dest) == 0
        @warn(
            "There are known correctness issues using Cbc with SOS " *
            "constraints and no binary variables.",
        )
    end
    return _index_map(src)
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.Utilities.UniversalFallback{OptimizerCache},
)
    attr = MOI.VariablePrimalStart()
    MOI.Utilities.throw_unsupported(
        src;
        excluded_attributes = Any[MOI.VariablePrimalStart()],
    )
    index_map = MOI.copy_to(dest, src.model)
    if attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        for (x_src, x_dest) in index_map.var_map
            value = MOI.get(src, attr, x_src)
            MOI.set(dest, attr, x_dest, value)
        end
    end
    return index_map
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    cache = MOI.default_cache(dest, Float64)
    src_cache = MOI.copy_to(cache, src)
    cache_dest = MOI.copy_to(dest, cache)
    index_map = MOI.IndexMap()
    for (src_x, cache_x) in src_cache.var_map
        index_map[src_x] = cache_dest[cache_x]
    end
    for (src_ci, cache_ci) in src_cache.con_map
        index_map[src_ci] = cache_dest[cache_ci]
    end
    return index_map
end

###
### supports and supports_constraint
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{
        <:Union{
            MOI.EqualTo{Float64},
            MOI.LessThan{Float64},
            MOI.GreaterThan{Float64},
            MOI.Interval{Float64},
            MOI.ZeroOne,
            MOI.Integer,
        },
    },
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.SOS1{Float64},MOI.SOS2{Float64}}},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{
        <:Union{
            MOI.EqualTo{Float64},
            MOI.LessThan{Float64},
            MOI.GreaterThan{Float64},
            MOI.Interval{Float64},
        },
    },
)
    return true
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end

###
### Variable starting values
###

function MOI.supports(
    ::Optimizer,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    x::MOI.VariableIndex,
    ::Nothing,
)
    delete!(model.variable_start, x)
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    x::MOI.VariableIndex,
    value,
)
    model.variable_start[x] = value
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    x::MOI.VariableIndex,
)
    return get(model.variable_start, x, nothing)
end

###
### Optimize and post-optimize functions
###

function MOI.optimize!(model::Optimizer)
    if !isempty(model.variable_start)
        columns = fill(Cint(0), length(model.variable_start))
        values = fill(0.0, length(model.variable_start))
        for (i, (variable, value)) in enumerate(model.variable_start)
            columns[i] = Cint(variable.value - 1)
            values[i] = value
        end
        Cbc_setMIPStartI(model, length(columns), columns, values)
    end
    t = time()
    model.variable_primal = nothing
    model.constraint_primal = nothing
    model.termination_status = Cbc_solve(model)
    model.has_solution = _result_count(model)
    model.solve_time = time() - t
    return
end

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return Cbc_getNumCols(model)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    if Cbc_getNumIntegers(model) == 0
        return MOI.get(model, MOI.ObjectiveValue())
    end
    return Cbc_getBestPossibleObjValue(model) + model.objective_constant
end

function MOI.get(model::Optimizer, ::MOI.NodeCount)::Int64
    return Cbc_getNodeCount(model)
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return Cbc_getObjValue(model) + model.objective_constant
end

# Cbc does not provide a native way of accessing the relative gap,
# use the Gurobi convention instead.
function MOI.get(model::Optimizer, ::MOI.RelativeGap)
    incumbent = MOI.get(model, MOI.ObjectiveValue())
    bound = MOI.get(model, MOI.ObjectiveBound())
    gap = abs(bound - incumbent) / abs(incumbent)
    return isnan(gap) ? Inf : gap
end

_update_cache(::Optimizer, data::Vector{Float64}, ::Any, ::Any) = data

function _update_cache(model::Optimizer, ::Nothing, f_p::F, f_n::G) where {F,G}
    p = f_p(model)
    n = f_n(model)
    if p == C_NULL
        return fill(NaN, n)
    end
    return unsafe_wrap(Array, p, (n,))
end

_get_cached_solution(data::Vector{Float64}, x) = data[x.value]

function _get_cached_solution(data::Vector{Float64}, x::Vector)
    return [data[xi.value] for xi in x]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    x::Union{MOI.VariableIndex,Vector{MOI.VariableIndex}},
)
    MOI.check_result_index_bounds(model, attr)
    model.variable_primal = _update_cache(
        model,
        model.variable_primal,
        Cbc_getColSolution,
        Cbc_getNumCols,
    )
    return _get_cached_solution(model.variable_primal, x)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
)
    MOI.check_result_index_bounds(model, attr)
    model.constraint_primal = _update_cache(
        model,
        model.constraint_primal,
        Cbc_getRowActivity,
        Cbc_getNumRows,
    )
    return _get_cached_solution(model.constraint_primal, c)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    index::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    MOI.check_result_index_bounds(model, attr)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(index.value))
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    sense = Cbc_getObjSense(model)
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

MOI.get(model::Optimizer, ::Status) = Cbc_status(model)
MOI.get(model::Optimizer, ::SecondaryStatus) = Cbc_secondaryStatus(model)

const _STATUS = Dict{Cint,String}(
    Cint(-1) => "before branchAndBound",
    Cint(0) =>
        "finished - check isProvenOptimal or isProvenInfeasible to see if solution found (or check value of best solution)",
    Cint(1) => "stopped - on maxnodes, maxsols, maxtime",
    Cint(2) => "execution abandoned due to numerical dificulties",
    Cint(5) => "user programmed interruption",
)

const _SECONDARY_STATUS = Dict{Cint,String}(
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
    Cbc_status          = $(_STATUS[model.termination_status])
    Cbc_secondaryStatus = $(_SECONDARY_STATUS[Cbc_secondaryStatus(model)])
    """
end

MOI.get(model::Optimizer, ::MOI.ResultCount) = model.has_solution ? 1 : 0

function _result_count(model::Optimizer)
    if Cbc_getNumIntegers(model) == 0
        # Cbc forwards the solve to the LP solver if there are no integers, so
        # check the termination status for the result count.
        return model.termination_status == 0 ? 1 : 0
    end
    return Cbc_numberSavedSolutions(model) > 0 ? 1 : 0
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    status = model.termination_status
    if status == -1
        return MOI.OPTIMIZE_NOT_CALLED
    elseif Cbc_isProvenOptimal(model) != 0
        return MOI.OPTIMAL
    elseif Cbc_isProvenInfeasible(model) != 0
        if Cbc_getNumIntegers(model) == 0
            # Why Cbc. For LPs, this could mean dual infeasible.
            return MOI.INFEASIBLE_OR_UNBOUNDED
        else
            return MOI.INFEASIBLE
        end
    elseif Cbc_isContinuousUnbounded(model) != 0
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif Cbc_isNodeLimitReached(model) != 0
        return MOI.NODE_LIMIT
    elseif Cbc_isSecondsLimitReached(model) != 0
        return MOI.TIME_LIMIT
    elseif Cbc_isSolutionLimitReached(model) != 0
        return MOI.SOLUTION_LIMIT
    elseif status == 0
        return MOI.OTHER_ERROR
    elseif status == 1
        return MOI.OTHER_LIMIT
    elseif status == 2
        return MOI.NUMERICAL_ERROR
    else
        @assert status == 5
        return MOI.INTERRUPTED
    end
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif MOI.get(model, MOI.ResultCount()) == 1
        return MOI.FEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(::Optimizer, ::MOI.DualStatus)
    return MOI.NO_SOLUTION
end
