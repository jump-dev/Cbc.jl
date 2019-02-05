using Compat.Test, MathOptInterface

const MOI  = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

@MOIU.model ModelForCachingOptimizer (MOI.ZeroOne, MOI.Integer) (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval) () () (MOI.SingleVariable,) (MOI.ScalarAffineFunction,) () ()

@testset "Continuous linear problems" begin
    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), Cbc.Optimizer())
    config = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals = false, infeas_certificates = false)
    MOIT.contlineartest(optimizer, config, [
    "linear1",  ## asks for ConstraintPrimal
    "linear2",  ## asks for ConstraintPrimal
    "linear7",  ## uses vector of constraints
    "linear10", ## asks for ConstraintPrimal
    "linear13", ## asks for ConstraintPrimal
    "linear14", ## asks for ConstraintPrimal
    "linear11", ## asks for ConstraintPrimal
    "linear15"  ## uses vector of constraints
    ])
end

@testset "Integer linear tests" begin

    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), Cbc.Optimizer())
    config = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals = false, infeas_certificates = false)
    # int1 excluded because asks for ConstraintPrimal
    # int2 excluded because uses vector of constraints
    MOIT.intlineartest(optimizer, config, ["int1", "int2"])

end


@testset "ModelLike tests" begin
    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), Cbc.Optimizer())

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
    optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), Cbc.Optimizer())

    MOIT.basic_constraint_tests(optimizer, config)

    MOIT.unittest(optimizer, config, [
        "solve_with_lowerbound", ## cannot get with strings
        "solve_blank_obj",
        "solve_qcp_edge_cases",
        "solve_qp_edge_cases",
        "solve_affine_interval", ## cannot get with strings
        "solve_affine_greaterthan", ## cannot get with strings
        "solve_with_upperbound",  ## cannot get with strings
        "solve_singlevariable_obj", ## cannot get with strings
        "solve_affine_equalto",  ## cannot get with strings
        "solve_affine_lessthan",  ## cannot get with strings
        "solve_constant_obj",  ## cannot get with strings
        "solve_affine_deletion_edge_cases", ## do not support vector of constraints
        ## TODO: fix new tests of objective edge cases
        "solve_duplicate_terms_obj",
        "solve_objbound_edge_cases"
    ])
end

@testset "Test params" begin
    default_model = ModelForCachingOptimizer{Float64}()
    x = MOI.add_variables(default_model, 100)
    terms = [MOI.ScalarAffineTerm(rand(), x[i]) for i in 1:length(x)]
    MOI.set(
        default_model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(terms, 0.0))
    MOI.set(default_model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    for i in 1:length(x)
        MOI.add_constraint(default_model, MOI.SingleVariable(x[i]), MOI.Integer())
        MOI.add_constraint(default_model, MOI.SingleVariable(x[i]), MOI.GreaterThan(0.0))
    end
    terms = [MOI.ScalarAffineTerm(rand(), x[i]) for i in 1:length(x)]
    MOI.add_constraint(
        default_model, MOI.ScalarAffineFunction(terms, 0.0), MOI.LessThan(31.0))

    model = Cbc.Optimizer(seconds=0)
    MOI.copy_to(model, default_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.TIME_LIMIT
    MOI.empty!(model)
    MOI.copy_to(model, default_model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.TIME_LIMIT
end
