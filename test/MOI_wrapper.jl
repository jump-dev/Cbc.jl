using Test, MathOptInterface

const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

import Cbc
const OPTIMIZER = Cbc.Optimizer()
MOI.set(OPTIMIZER, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "COIN Branch-and-Cut (Cbc)"
end

@testset "supports_default_copy_to" begin
    @test !MOIU.supports_allocate_load(OPTIMIZER, false)
    @test !MOIU.supports_allocate_load(OPTIMIZER, true)
    @test !MOIU.supports_default_copy_to(OPTIMIZER, false)
    @test !MOIU.supports_default_copy_to(OPTIMIZER, true)
end

const CACHE = MOIU.UniversalFallback(MOIU.Model{Float64}())
const CACHED = MOIU.CachingOptimizer(CACHE, OPTIMIZER)
const BRIDGED = MOI.Bridges.full_bridge_optimizer(CACHED, Float64)

const CONFIG = MOIT.TestConfig(duals = false, infeas_certificates = false)

@testset "basic_constraint_tests" begin
    MOIT.basic_constraint_tests(CACHED, CONFIG)
end

@testset "Unit" begin
    MOIT.unittest(BRIDGED, CONFIG, [
        "number_threads", # FIXME implement `MOI.NumberOfThreads`
        "solve_time", # FIXME implement `MOI.SolveTime`
        "solve_unbounded_model", # INFEASIBLE_OR_UNBOUNDED instead of DUAL_INFEASIBLE
        "delete_soc_variables", "solve_qcp_edge_cases", "solve_qp_edge_cases"  # No quadratics
    ])
end

@testset "ModelLike" begin
    @test MOI.get(CACHED, MOI.SolverName()) == "COIN Branch-and-Cut (Cbc)"
    @testset "default_objective_test" begin
         MOIT.default_objective_test(CACHED)
     end
     @testset "default_status_test" begin
         MOIT.default_status_test(CACHED)
     end
    @testset "nametest" begin
        MOIT.nametest(CACHED)
    end
    @testset "validtest" begin
        MOIT.validtest(CACHED)
    end
    @testset "emptytest" begin
        # Requires VectorOfVariables
        # MOIT.emptytest(CACHED)
    end
    @testset "orderedindicestest" begin
        MOIT.orderedindicestest(CACHED)
    end
    @testset "copytest" begin
        # Requires VectorOfVariables
        # MOIT.copytest(CACHED, MOIU.CachingOptimizer(
        #     ModelForCachingOptimizer{Float64}(),
        #     Cbc.Optimizer()
        # ))
    end
end

@testset "Continuous Linear" begin
    MOIT.contlineartest(BRIDGED, CONFIG)
end

@testset "Integer Linear" begin
    MOIT.intlineartest(BRIDGED, CONFIG, [
        "indicator1", "indicator2", "indicator3"
    ])
end

@testset "Test params" begin
    knapsack_model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(knapsack_model, """
        variables: x, y
        maxobjective: x + y
        c1: x in ZeroOne()
        c2: y in ZeroOne()
        c3: x + y <= 1.0
    """)
    model = Cbc.Optimizer(maxNodes = 0, presolve = "off", cuts = "off",
                          heur = "off", logLevel = 0)
    MOI.copy_to(model, knapsack_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.NODE_LIMIT
    # We also check that options are not destroyed on `empty!`.
    MOI.empty!(model)
    MOI.copy_to(model, knapsack_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.NODE_LIMIT
end
