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
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            Cbc.Optimizer(),
        ),
        Float64,
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
            # TODO(odow): bug in Cbc.jl
            "test_constraint_Indicator_ACTIVATE_ON_ZERO",
            "test_model_copy_to_UnsupportedAttribute",
            "test_model_ModelFilter_AbstractConstraintAttribute",
            "test_objective_FEASIBILITY_SENSE_clears_objective",
            # TODO(odow): bug in MOI
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
            # TODO(odow): upstream bug in Cbc
            "test_linear_Indicator_",
            "test_linear_SOS1_integration",
            "test_linear_SOS2_integration",
            "test_solve_SOS2_",
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
    model = Cbc.Optimizer(
        maxSol = 1,
        presolve = "off",
        cuts = "off",
        heur = "off",
        logLevel = 0,
    )
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

end

TestMOIWrapper.runtests()
