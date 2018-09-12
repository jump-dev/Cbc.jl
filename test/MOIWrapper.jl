using Compat.Test, MathOptInterface, MathOptInterface.Test, MathOptInterface.Utilities

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test
const MOIU = MathOptInterface.Utilities
const MOIB = MathOptInterface.Bridges

@MOIU.model ModelForCachingOptimizer (MOI.ZeroOne, MOI.Integer) (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval) () () (MOI.SingleVariable,) (MOI.ScalarAffineFunction,) () ()

@testset "Continuous linear problems" begin
    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), CbcOptimizer())
    config = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals = false, infeas_certificates = false)
    MOIT.contlineartest(optimizer, config, [
    "linear1",  ## asks for ConstraintPrimal
    "linear2",  ## asks for ConstraintPrimal
    "linear7",  ## uses vector of constraints
    "linear10", ## asks for ConstraintPrimal
    "linear13", ## asks for ConstraintPrimal
    "linear14", ## asks for ConstraintPrimal
    "linear11", ## fails due to a bug in Clp.
    "linear15"  ## uses vector of constraints
    ])
end

@testset "Integer linear tests" begin

    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), CbcOptimizer())
    config = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals = false, infeas_certificates = false)
    # int1 excluded because asks for ConstraintPrimal
    # int2 excluded because uses vector of constraints
    MOIT.intlineartest(optimizer, config, ["int1", "int2"])

end


@testset "ModelLike tests" begin
    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), CbcOptimizer())

    # MOIT.nametest(optimizer) ## uses names/strings
    # @testset "validtest" begin ## at some moment inside the test it asks for the number of constraints of a type that is not supported by my model and caching optimizer returns an error. I cannot return the number of constraints of a specific type
    #     MOIT.validtest(optimizer)
    # end
    # @testset "emptytest" begin ## vector of constraints
    #     MOIT.emptytest(optimizer)
    # end
    @testset "orderedindicestest" begin
        MOIT.orderedindicestest(optimizer)
    end
    # @testset "canaddconstrainttest" begin ## do not pass due to vector of variables
    #     MOIT.canaddconstrainttest(optimizer, Float64, Complex{Float64})
    # end
    # @testset "copytest" begin ## do not pass due to vector of variables
    #     MOIT.copytest(optimizer.optimizer, optimizer)
    # end
end

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), CbcOptimizer())

    MOIT.basic_constraint_tests(optimizer, config)

    MOIT.unittest(optimizer, config, [
        "solve_with_lowerbound", ## cannot get with strings
        "solve_blank_obj",
        "solve_qcp_edge_cases",
        "solve_qp_edge_cases",
        "solve_affine_interval", ## cannot get with strings
        "solve_affine_greaterthan", ## cannot get with strings
        "solve_with_upperbound",  ## cannot get with strings
        "getvariable",  ## cannot get with strings
        "solve_singlevariable_obj", ## cannot get with strings
        "solve_affine_equalto",  ## cannot get with strings
        "getconstraint", ## cannot get constraint
        "solve_affine_lessthan",  ## cannot get with strings
        "solve_constant_obj",  ## cannot get with strings
        "solve_affine_deletion_edge_cases", ## do not support vector of constraints
        ## TODO: fix new tests of objective edge cases
        "solve_duplicate_terms_obj", 
        "solve_objbound_edge_cases"
    ])
end
