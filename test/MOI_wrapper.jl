using Cbc, Test, MathOptInterface

const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

@MOIU.model(ModelForCachingOptimizer,
    (MOI.ZeroOne, MOI.Integer),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval),
    (),
    (),
    (MOI.SingleVariable,),
    (MOI.ScalarAffineFunction,),
    (),
    ()
)

const OPTIMIZER = MOIU.CachingOptimizer(
    MOIU.UniversalFallback(ModelForCachingOptimizer{Float64}()),
    Cbc.Optimizer(logLevel = 0)
)

const CONFIG = MOIT.TestConfig(duals = false, infeas_certificates = false)

@testset "basic_constraint_tests" begin
    MOIT.basic_constraint_tests(OPTIMIZER, CONFIG)
end

@testset "Unit Tests" begin
    MOIT.unittest(OPTIMIZER, CONFIG, [
        "solve_affine_deletion_edge_cases",  # VectorAffineFunction
        "solve_duplicate_terms_vector_affine",  # VectorAffineFunction
        "solve_qcp_edge_cases", "solve_qp_edge_cases"  # No quadratics
    ])
end

@testset "ModelLike tests" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "COIN Branch-and-Cut (Cbc)"
    @testset "default_objective_test" begin
         MOIT.default_objective_test(OPTIMIZER)
     end
     @testset "default_status_test" begin
         MOIT.default_status_test(OPTIMIZER)
     end
    @testset "nametest" begin
        MOIT.nametest(OPTIMIZER)
    end
    @testset "validtest" begin
        MOIT.validtest(OPTIMIZER)
    end
    @testset "emptytest" begin
        # Requires VectorOfVariables
        # MOIT.emptytest(OPTIMIZER)
    end
    @testset "orderedindicestest" begin
        MOIT.orderedindicestest(OPTIMIZER)
    end
    @testset "copytest" begin
        # Requires VectorOfVariables
        # MOIT.copytest(OPTIMIZER, MOIU.CachingOptimizer(
        #     ModelForCachingOptimizer{Float64}(),
        #     Cbc.Optimizer()
        # ))
    end
end

@testset "contlineartest" begin
    MOIT.contlineartest(OPTIMIZER, CONFIG, [
        "linear7", "linear15"  # VectorAffineFunction
    ])
end

@testset "intlineartest" begin
    MOIT.intlineartest(OPTIMIZER, CONFIG, [
        "int2"  # Requires Special-Ordered-Sets
    ])
end

@testset "Test params" begin
    knapsack_model = ModelForCachingOptimizer{Float64}()
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
