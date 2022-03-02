module TestMOIWrapper

using Test
using MathOptInterface
import Cbc

const MOI = MathOptInterface

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_SolverName()
    @test MOI.get(Cbc.Optimizer(), MOI.SolverName()) ==
          "COIN Branch-and-Cut (Cbc)"
end

function test_supports_incremental_interface()
    @test !MOI.supports_incremental_interface(Cbc.Optimizer())
    return
end

function test_runtests()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(Cbc.Optimizer; with_bridge_type = Float64),
    )
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            exclude = Any[
                MOI.ConstraintDual,
                MOI.DualObjectiveValue,
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
            ],
        ),
        exclude = [
            # TODO(odow): upstream bug in Cbc
            "test_linear_Indicator_",
            "test_linear_SOS1_integration",
            "test_linear_SOS2_integration",
            "test_solve_SOS2_add_and_delete",
            # Can't prove infeasible.
            "test_conic_NormInfinityCone_INFEASIBLE",
            "test_conic_NormOneCone_INFEASIBLE",
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
        ],
    )
    return
end

function test_params()
    # Note: we generate a non-trivial problem to ensure that Cbc struggles to
    # find a solution at the root node.
    knapsack_model = MOI.Utilities.Model{Float64}()
    N = 100
    x = MOI.add_variables(knapsack_model, N)
    MOI.add_constraint.(knapsack_model, x, MOI.ZeroOne())
    MOI.add_constraint(
        knapsack_model,
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([1 + sin(i) for i in 1:N], x),
            0.0,
        ),
        MOI.LessThan(10.0),
    )
    MOI.set(
        knapsack_model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([cos(i) for i in 1:N], x),
            0.0,
        ),
    )
    model = Cbc.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute("maxSol"), 1)
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("cuts"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("heur"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("logLevel"), 0)
    MOI.copy_to(model, knapsack_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.SOLUTION_LIMIT
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    MOI.empty!(model)
    MOI.copy_to(model, knapsack_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.SOLUTION_LIMIT
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    return
end

"""
    test_threads()

Test solving a model with the threads parameter set.

See issues #112 and #186.
"""
function test_threads()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(Cbc.Optimizer; with_bridge_type = Float64),
    )
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("cuts"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("heur"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("threads"), 4)
    MOI.set(model, MOI.RawOptimizerAttribute("logLevel"), 3)
    N = 100
    x = MOI.add_variables(model, N)
    MOI.add_constraint.(model, x, MOI.ZeroOne())
    w = [1 + sin(i) for i in 1:N]
    c = [1 + cos(i) for i in 1:N]
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(w, x), 0.0),
        MOI.LessThan(10.0),
    )
    obj = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0)
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    return
end

function test_PrimalStatus()
    model = MOI.Utilities.Model{Float64}()
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.GreaterThan(1.0))
    MOI.add_constraint(model, x, MOI.LessThan(0.0))
    cbc = Cbc.Optimizer()
    MOI.copy_to(cbc, model)
    MOI.optimize!(cbc)
    MOI.get(cbc, MOI.PrimalStatus()) == MOI.NO_SOLUTION
    return
end

# The test_linear_SOS1_integration test with the additional requirement that all
# variables are integer.
function test_SOS1()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            Cbc.Optimizer(),
        ),
        Float64,
    )
    config = MOI.Test.Config()
    MOI.set(model, MOI.Silent(), true)
    @test MOI.supports_constraint(
        model,
        MOI.VectorOfVariables,
        MOI.SOS1{Float64},
    )
    @test MOI.supports_constraint(
        model,
        MOI.VariableIndex,
        MOI.LessThan{Float64},
    )
    v = MOI.add_variables(model, 3)
    MOI.add_constraint.(model, v, MOI.Integer())
    @test MOI.get(model, MOI.NumberOfVariables()) == 3
    vc1 = MOI.add_constraint(model, v[1], MOI.LessThan(1.0))
    @test vc1.value == v[1].value
    vc2 = MOI.add_constraint(model, v[2], MOI.LessThan(1.0))
    @test vc2.value == v[2].value
    vc3 = MOI.add_constraint(model, v[3], MOI.LessThan(2.0))
    @test vc3.value == v[3].value
    c1 = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([v[1], v[2]]),
        MOI.SOS1([1.0, 2.0]),
    )
    c2 = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([v[1], v[3]]),
        MOI.SOS1([1.0, 2.0]),
    )
    @test MOI.get(
        model,
        MOI.NumberOfConstraints{MOI.VectorOfVariables,MOI.SOS1{Float64}}(),
    ) == 2
    #=
        To allow for permutations in the sets and variable vectors
        we're going to sort according to the weights
    =#
    cs_sos = MOI.get(model, MOI.ConstraintSet(), c2)
    cf_sos = MOI.get(model, MOI.ConstraintFunction(), c2)
    p = sortperm(cs_sos.weights)
    @test isapprox(cs_sos.weights[p], [1.0, 2.0], config)
    @test cf_sos.variables[p] == v[[1, 3]]
    objf =
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0, 1.0, 1.0], v), 0.0)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        objf,
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    @test MOI.get(model, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
    @test MOI.get(model, MOI.ResultCount()) >= 1
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 3, config)
    @test isapprox(MOI.get(model, MOI.VariablePrimal(), v), [0, 1, 2], config)
    MOI.delete(model, c1)
    MOI.delete(model, c2)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
    @test MOI.get(model, MOI.ResultCount()) >= 1
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 5, config)
    @test isapprox(MOI.get(model, MOI.VariablePrimal(), v), [1, 1, 2], config)
    return
end

# The test_linear_SOS2_integration test with the additional requirement that all
# variables are integer.
function test_SOS2()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            Cbc.Optimizer(),
        ),
        Float64,
    )
    config = MOI.Test.Config()
    MOI.set(model, MOI.Silent(), true)
    @test MOI.supports_constraint(
        model,
        MOI.VectorOfVariables,
        MOI.SOS1{Float64},
    )
    @test MOI.supports_constraint(
        model,
        MOI.VectorOfVariables,
        MOI.SOS2{Float64},
    )
    v = MOI.add_variables(model, 10)
    @test MOI.get(model, MOI.NumberOfVariables()) == 10
    bin_constraints = []
    for i in 1:8
        vc = MOI.add_constraint(model, v[i], MOI.Interval(0.0, 2.0))
        @test vc.value == v[i].value
        push!(bin_constraints, MOI.add_constraint(model, v[i], MOI.ZeroOne()))
        @test bin_constraints[i].value == v[i].value
    end
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([1.0, 2.0, 3.0, -1.0], v[[1, 2, 3, 9]]),
            0.0,
        ),
        MOI.EqualTo(0.0),
    )
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.(
                [5.0, 4.0, 7.0, 2.0, 1.0, -1.0],
                v[[4, 5, 6, 7, 8, 10]],
            ),
            0.0,
        ),
        MOI.EqualTo(0.0),
    )
    MOI.add_constraint(
        model,
        MOI.VectorOfVariables(v[[1, 2, 3]]),
        MOI.SOS1([1.0, 2.0, 3.0]),
    )
    vv = MOI.VectorOfVariables(v[[4, 5, 6, 7, 8]])
    sos2 = MOI.SOS2([5.0, 4.0, 7.0, 2.0, 1.0])
    c = MOI.add_constraint(model, vv, sos2)
    #=
        To allow for permutations in the sets and variable vectors
        we're going to sort according to the weights
    =#
    cs_sos = MOI.get(model, MOI.ConstraintSet(), c)
    cf_sos = MOI.get(model, MOI.ConstraintFunction(), c)
    p = sortperm(cs_sos.weights)
    @test isapprox(cs_sos.weights[p], [1.0, 2.0, 4.0, 5.0, 7.0], config)
    @test cf_sos.variables[p] == v[[8, 7, 5, 4, 6]]
    objf = MOI.ScalarAffineFunction(
        MOI.ScalarAffineTerm.([1.0, 1.0], [v[9], v[10]]),
        0.0,
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        objf,
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    @test MOI.get(model, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
    @test MOI.get(model, MOI.ResultCount()) >= 1
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 15.0, config)
    @test isapprox(
        MOI.get(model, MOI.VariablePrimal(), v),
        [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 3.0, 12.0],
        config,
    )
    for cref in bin_constraints
        MOI.delete(model, cref)
    end
    MOI.add_constraint.(model, v, MOI.Integer())
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
    @test MOI.get(model, MOI.ResultCount()) >= 1
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 30.0, config)
    @test isapprox(
        MOI.get(model, MOI.VariablePrimal(), v),
        [0.0, 0.0, 2.0, 2.0, 0.0, 2.0, 0.0, 0.0, 6.0, 24.0],
        config,
    )
    return
end

"""
    test_VariablePrimalStart()

Testing that VariablePrimalStart is actually applied is a little convoluted.

We formulate a MIP with various setttings turned off to avoid a trivial solve in
presolve.

Then we solve and return the optimal primal solution and the number of nodes
visited.

For the second pass, we rebuild the same MIP, but this time we pass the optimal
solution as the VariablePrimalStart, and we set maxSol=1 to force Cbc to exit
after finding a single solution. Because we passed a primal feasible point, it
should return the optimal solution after exploring 0 nodes.
"""
function test_VariablePrimalStart()
    function formulate_and_solve(start)
        model = MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            MOI.instantiate(Cbc.Optimizer; with_bridge_type = Float64),
        )
        MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
        MOI.set(model, MOI.RawOptimizerAttribute("cuts"), "off")
        MOI.set(model, MOI.RawOptimizerAttribute("heur"), "off")
        MOI.set(model, MOI.RawOptimizerAttribute("logLevel"), 0)
        N = 100
        x = MOI.add_variables(model, N)
        MOI.add_constraint.(model, x, MOI.ZeroOne())
        w = [1 + sin(i) for i in 1:N]
        c = [1 + cos(i) for i in 1:N]
        if start !== nothing
            MOI.set(model, MOI.RawOptimizerAttribute("maxSol"), 1)
            MOI.set.(model, MOI.VariablePrimalStart(), x, start)
        end
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(w, x), 0.0),
            MOI.LessThan(10.0),
        )
        obj = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0)
        MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
        MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.optimize!(model)
        sol = MOI.get.(model, MOI.VariablePrimal(), x)
        return sol, MOI.get(model, MOI.NodeCount())
    end
    x_sol, nodes = formulate_and_solve(nothing)
    y_sol, nodes_start = formulate_and_solve(x_sol)
    @test x_sol == y_sol
    @test nodes > 0
    @test nodes_start == 0
    return
end

end

TestMOIWrapper.runtests()
