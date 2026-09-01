# Copyright (c) 2013: Cbc.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestMOIWrapper

using Test

import Cbc
import MathOptInterface as MOI

function runtests()
    is_test(name) = startswith("$(name)", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_SolverName()
    @test MOI.get(Cbc.Optimizer(), MOI.SolverName()) ==
          "COIN Branch-and-Cut (Cbc)"
    return
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
    @test MOI.get(model, MOI.RawOptimizerAttribute("maxSol")) == "1"
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
    @test MOI.get(model, MOI.RelativeGap()) >= 0
    @test MOI.is_set_by_optimize(Cbc.Status())
    @test MOI.get(model, Cbc.Status()) == 1
    @test MOI.is_set_by_optimize(Cbc.SecondaryStatus())
    @test MOI.get(model, Cbc.SecondaryStatus()) == 6
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

function test_issue_187()
    if true
        # This test segfaults (unreliably) on all platforms
        @test_broken 1 == 2
        return
    end
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(Cbc.Optimizer; with_bridge_type = Float64),
    )
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint.(model, x, MOI.ZeroOne())
    @test MOI.get.(model, MOI.VariablePrimalStart(), x) == [nothing, nothing]
    MOI.set.(model, MOI.VariablePrimalStart(), x, 0.0)
    @test MOI.get.(model, MOI.VariablePrimalStart(), x) == [0.0, 0.0]
    y = MOI.add_variables(model, 2)
    MOI.add_constraint.(model, y, MOI.ZeroOne())

    MOI.add_constraint(
        model,
        MOI.Utilities.operate(vcat, Float64, x[1], 1.0 * y[1]),
        MOI.Indicator{MOI.ACTIVATE_ON_ONE}(MOI.EqualTo(1.0)),
    )
    MOI.add_constraint(
        model,
        MOI.Utilities.operate(vcat, Float64, x[2], 1.0 * y[1]),
        MOI.Indicator{MOI.ACTIVATE_ON_ONE}(MOI.EqualTo(0.0)),
    )
    MOI.set.(model, MOI.VariablePrimalStart(), y, 0.0)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0)
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test ≈(sum(MOI.get(model, MOI.VariablePrimal(), x)), 1.0, atol = 1e-4)
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

function test_variable_name()
    for (name, inner) in [("abc", "abc"), ("απ", "C0000000")]
        model = MOI.Utilities.Model{Float64}()
        x = MOI.add_variable(model)
        MOI.set(model, MOI.VariableName(), x, name)
        cbc = Cbc.Optimizer()
        @test !MOI.supports(cbc, MOI.VariableName(), MOI.VariableIndex)
        MOI.set(cbc, Cbc.SetVariableNames(), true)
        @test MOI.supports(cbc, MOI.VariableName(), MOI.VariableIndex)
        index_map = MOI.copy_to(cbc, model)
        @test MOI.get(cbc, MOI.VariableName(), index_map[x]) == inner
    end
    return
end

"""
This example segfaults if variable names are set. Turn them off and check we get
a feasible solution.
"""
function test_segfault()
    src = MOI.FileFormats.MOF.Model()
    MOI.read_from_file(src, joinpath(@__DIR__, "segfault.mof.json"))
    cbc = Cbc.Optimizer()
    @test MOI.supports(cbc, Cbc.SetVariableNames())
    @test MOI.get(cbc, Cbc.SetVariableNames()) == false
    MOI.set(cbc, Cbc.SetVariableNames(), true)
    @test MOI.get(cbc, Cbc.SetVariableNames()) == true
    MOI.set(cbc, Cbc.SetVariableNames(), false)
    @test MOI.get(cbc, Cbc.SetVariableNames()) == false
    index_map = MOI.copy_to(cbc, src)
    MOI.optimize!(cbc)
    @test MOI.get(cbc, MOI.TerminationStatus()) == MOI.OPTIMAL
    return
end

function test_get_objective_sense()
    for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE, MOI.FEASIBILITY_SENSE)
        model = Cbc.Optimizer()
        src = MOI.Utilities.Model{Float64}()
        MOI.set(src, MOI.ObjectiveSense(), sense)
        MOI.copy_to(model, src)
        @test MOI.get(model, MOI.ObjectiveSense()) == sense
    end
    return
end

end

TestMOIWrapper.runtests()
