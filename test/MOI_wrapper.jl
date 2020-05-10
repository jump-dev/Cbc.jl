using Test
using MathOptInterface
import Cbc

const MOI = MathOptInterface

const OPTIMIZER = Cbc.Optimizer()
MOI.set(OPTIMIZER, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "COIN Branch-and-Cut (Cbc)"
end

@testset "supports_default_copy_to" begin
    @test !MOI.Utilities.supports_allocate_load(OPTIMIZER, false)
    @test !MOI.Utilities.supports_allocate_load(OPTIMIZER, true)
    @test !MOI.Utilities.supports_default_copy_to(OPTIMIZER, false)
    @test !MOI.Utilities.supports_default_copy_to(OPTIMIZER, true)
end

const CACHED = MOI.Utilities.CachingOptimizer(
    MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    OPTIMIZER
)
const BRIDGED = MOI.Bridges.full_bridge_optimizer(CACHED, Float64)

const CONFIG = MOI.Test.TestConfig(duals = false, infeas_certificates = false)

@testset "basic_constraint_tests" begin
    MOI.Test.basic_constraint_tests(CACHED, CONFIG)
end

@testset "Unit" begin
    MOI.Test.unittest(BRIDGED, CONFIG, [
        # TODO(odow): implement attributes:
        "number_threads",

        # INFEASIBLE_OR_UNBOUNDED instead of DUAL_INFEASIBLE
        "solve_unbounded_model",

        # No quadratics
        "delete_soc_variables",
        "solve_qcp_edge_cases",
        "solve_qp_edge_cases",
    ])
end

@testset "ModelLike" begin
    @test MOI.get(CACHED, MOI.SolverName()) == "COIN Branch-and-Cut (Cbc)"
    @testset "default_objective_test" begin
         MOI.Test.default_objective_test(CACHED)
     end
     @testset "default_status_test" begin
         MOI.Test.default_status_test(CACHED)
     end
    @testset "nametest" begin
        MOI.Test.nametest(CACHED)
    end
    @testset "validtest" begin
        MOI.Test.validtest(CACHED)
    end
    @testset "emptytest" begin
        MOI.Test.emptytest(BRIDGED)
    end
    @testset "orderedindicestest" begin
        MOI.Test.orderedindicestest(CACHED)
    end
end

@testset "Continuous Linear" begin
    MOI.Test.contlineartest(BRIDGED, CONFIG)
end

@testset "Integer Linear" begin
    MOI.Test.intlineartest(BRIDGED, CONFIG, [
        # Cbc does not support indicator constraints.
        "indicator1", "indicator2", "indicator3", "indicator4",
        # TODO(odow): needs MOI at least 0.9.14.
        "semiconttest", "semiinttest"
    ])
end

@testset "Test params" begin
    # Note: we generate a non-trivial problem to ensure that Cbc struggles to
    # find a solution at the root node.
    knapsack_model = MOI.Utilities.Model{Float64}()
    N = 100
    x = MOI.add_variables(knapsack_model, N)
    MOI.add_constraint.(knapsack_model, MOI.SingleVariable.(x), MOI.ZeroOne())
    MOI.add_constraint(
        knapsack_model,
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([1 + sin(i) for i = 1:N], x),
            0.0
        ),
        MOI.LessThan(10.0),
    )
    MOI.set(
        knapsack_model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([cos(i) for i = 1:N], x),
            0.0
        )
    )
    model = Cbc.Optimizer(
        maxNodes = 0, presolve = "off", cuts = "off", heur = "off", logLevel = 0
    )
    MOI.copy_to(model, knapsack_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.NODE_LIMIT
    MOI.empty!(model)
    MOI.copy_to(model, knapsack_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.NODE_LIMIT
end
