# using MathOptInterface
# const MOI = MathOptInterface
# const MOIT = MOI.Test
# const MOIB = MOI.Bridges
#
# const MOIU = MOI.Utilities
# MOIU.@model SCSModelData () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, ExponentialCone, PositiveSemidefiniteConeTriangle) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)
# const optimizer = MOIU.CachingOptimizer(SCSModelData{Float64}(), SCSOptimizer())
#
# # linear9test needs 1e-3 with SCS < 2.0 and 5e-1 with SCS 2.0
# # linear2test needs 1e-4
# const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)
#
# @testset "Continuous linear problems" begin
#     # AlmostSuccess for linear9 with SCS 2
#     MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config, ["linear9"])
# end
#








using Base.Test, MathOptInterface, MathOptInterface.Test, MathOptInterface.Utilities

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test
const MOIU = MathOptInterface.Utilities

# @testset "Cont Linear tests" begin
#     linconfig = MOIT.TestConfig(modify_lhs = false, duals = false)
#     @testset "Default Solver"  begin
#         solver = CbcOptimizer(LogLevel = 0)
#         MOIT.contlineartest(solver, linconfig, ["linear10","linear12","linear8a","linear8b","linear8c"])
#     end
# end

@MOIU.model ModelForCachingOptimizer (ZeroOne, Integer) (EqualTo, GreaterThan, LessThan, Interval) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, RotatedSecondOrderCone, GeometricMeanCone, ExponentialCone, DualExponentialCone, PositiveSemidefiniteConeTriangle, RootDetConeTriangle, LogDetConeTriangle) () (SingleVariable,) (ScalarAffineFunction,ScalarQuadraticFunction) (VectorOfVariables,) (VectorAffineFunction,)

@testset "Single variable EqualTo" begin
    m = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), MOIU.Manual)

    v = MOI.addvariable!(m)

    @test MOI.canaddconstraint(m, MOI.SingleVariable, MOI.EqualTo{Float64})
    vc = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.EqualTo(5.0))

    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([3.0], [v]), 0.0)
    @test MOI.canset(m, MOI.ObjectiveFunction{typeof(saf)}())
    MOI.set!(m, MOI.ObjectiveFunction{typeof(saf)}(), saf)
    @test MOI.get(m, MOIU.AttributeFromModelCache(MOI.ObjectiveFunction{typeof(saf)}())) ≈ saf
    @test MOI.get(m, MOI.ObjectiveFunction{typeof(saf)}()) ≈ saf

    solver = CbcOptimizer()

    MOI.copy!(solver, m)
    myStatus = MOI.optimize!(solver)

    @test MOI.get(solver, MOI.ObjectiveValue()) == 15.0
    @test MOI.get(solver, MOI.ObjectiveBound()) == 15.0

end

@testset "Single variable GreaterThan" begin
    m = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), MOIU.Manual)

    v = MOI.addvariable!(m)

    @test MOI.canaddconstraint(m, MOI.SingleVariable, MOI.GreaterThan{Float64})
    vc2 = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.GreaterThan(5.0))

    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([3.0], [v]), 0.0)
    @test MOI.canset(m, MOI.ObjectiveFunction{typeof(saf)}())
    MOI.set!(m, MOI.ObjectiveFunction{typeof(saf)}(), saf)
    @test MOI.get(m, MOIU.AttributeFromModelCache(MOI.ObjectiveFunction{typeof(saf)}())) ≈ saf
    @test MOI.get(m, MOI.ObjectiveFunction{typeof(saf)}()) ≈ saf

    solver = CbcOptimizer()

    MOI.copy!(solver, m)
    myStatus = MOI.optimize!(solver)

    @test MOI.get(solver, MOI.ObjectiveValue()) == 15.0
    @test MOI.get(solver, MOI.ObjectiveBound()) == 15.0

end


@testset "Single variable LessThan" begin
    m = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), MOIU.Manual)

    v = MOI.addvariable!(m)

    @test MOI.canaddconstraint(m, MOI.SingleVariable, MOI.LessThan{Float64})
    vc2 = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.LessThan(5.0))

    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([3.0], [v]), 0.0)
    @test MOI.canset(m, MOI.ObjectiveFunction{typeof(saf)}())
    MOI.set!(m, MOI.ObjectiveFunction{typeof(saf)}(), saf)
    @test MOI.get(m, MOIU.AttributeFromModelCache(MOI.ObjectiveFunction{typeof(saf)}())) ≈ saf
    @test MOI.get(m, MOI.ObjectiveFunction{typeof(saf)}()) ≈ saf

    solver = CbcOptimizer()

    MOI.copy!(solver, m)
    MOI.set!(solver, MOI.ObjectiveSense(), MOI.MaxSense)
    @test MOI.get(solver, MOI.ObjectiveSense()) == MOI.MaxSense

    myStatus = MOI.optimize!(solver)

    @test MOI.get(solver, MOI.ObjectiveValue()) == 15.0
    @test MOI.get(solver, MOI.ObjectiveBound()) == 15.0

end


@testset "Single negative variable" begin
    m = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), MOIU.Manual)

    v = MOI.addvariable!(m)

    @test MOI.canaddconstraint(m, MOI.SingleVariable, MOI.LessThan{Float64})
    vc2 = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.GreaterThan(-1.0))

    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [v]), 0.0)
    @test MOI.canset(m, MOI.ObjectiveFunction{typeof(saf)}())
    MOI.set!(m, MOI.ObjectiveFunction{typeof(saf)}(), saf)
    @test MOI.get(m, MOIU.AttributeFromModelCache(MOI.ObjectiveFunction{typeof(saf)}())) ≈ saf
    @test MOI.get(m, MOI.ObjectiveFunction{typeof(saf)}()) ≈ saf

    solver = CbcOptimizer()

    MOI.copy!(solver, m)
    MOI.set!(solver, MOI.ObjectiveSense(), MOI.MinSense)
    @test MOI.get(solver, MOI.ObjectiveSense()) == MOI.MinSense

    myStatus = MOI.optimize!(solver)

    @test MOI.get(solver, MOI.ObjectiveValue()) == -1.0
    @test MOI.get(solver, MOI.ObjectiveBound()) == -1.0

end



@testset "Single Integer max sense" begin
    m = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), MOIU.Manual)

    v = MOI.addvariable!(m)

    @test MOI.canaddconstraint(m, MOI.SingleVariable, MOI.LessThan{Float64})
    vc1 = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.LessThan(1.5))
    vc2 = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.Integer())

    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([3.0], [v]), 0.0)
    @test MOI.canset(m, MOI.ObjectiveFunction{typeof(saf)}())
    MOI.set!(m, MOI.ObjectiveFunction{typeof(saf)}(), saf)
    @test MOI.get(m, MOIU.AttributeFromModelCache(MOI.ObjectiveFunction{typeof(saf)}())) ≈ saf
    @test MOI.get(m, MOI.ObjectiveFunction{typeof(saf)}()) ≈ saf

    solver = CbcOptimizer()

    MOI.copy!(solver, m)
    MOI.set!(solver, MOI.ObjectiveSense(), MOI.MaxSense)
    @test MOI.get(solver, MOI.ObjectiveSense()) == MOI.MaxSense
    myStatus = MOI.optimize!(solver)

    @test MOI.get(solver, MOI.ObjectiveValue()) == 3.0
    @test MOI.get(solver, MOI.ObjectiveBound()) == 3.0
    @test MOI.get(solver, MOI.VariablePrimal(), MOI.VariableIndex(1)) == 1

end


@testset "Single Integer min sense" begin
    m = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), MOIU.Manual)

    v = MOI.addvariable!(m)

    @test MOI.canaddconstraint(m, MOI.SingleVariable, MOI.GreaterThan{Float64})
    vc1 = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.GreaterThan(0.5))
    vc2 = MOI.addconstraint!(m, MOI.SingleVariable(v), MOI.Integer())

    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([3.0], [v]), 0.0)
    @test MOI.canset(m, MOI.ObjectiveFunction{typeof(saf)}())
    MOI.set!(m, MOI.ObjectiveFunction{typeof(saf)}(), saf)
    @test MOI.get(m, MOIU.AttributeFromModelCache(MOI.ObjectiveFunction{typeof(saf)}())) ≈ saf
    @test MOI.get(m, MOI.ObjectiveFunction{typeof(saf)}()) ≈ saf

    solver = CbcOptimizer()

    MOI.copy!(solver, m)
    MOI.set!(solver, MOI.ObjectiveSense(), MOI.MinSense)
    @test MOI.get(solver, MOI.ObjectiveSense()) == MOI.MinSense

    myStatus = MOI.optimize!(solver)

    @test MOI.get(solver, MOI.ObjectiveValue()) == 3.0
    @test MOI.get(solver, MOI.ObjectiveBound()) == 3.0
    @test MOI.get(solver, MOI.VariablePrimal(), MOI.VariableIndex(1)) == 1

end


@testset "scalar affine constraint" begin
    m = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), MOIU.Manual)

    vars1 = MOI.addvariables!(m, 3)


    MOI.addconstraint!(m, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0, 1.0], vars1), 0.0), MOI.LessThan(1.0))

    # MOI.addconstraint!(m, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0, 1.0], vars1), 0.0), MOI.GreaterThan(1.0))

    vc1 = MOI.addconstraint!(m, MOI.SingleVariable(vars1[1]), MOI.GreaterThan(0.0))
    vc2 = MOI.addconstraint!(m, MOI.SingleVariable(vars1[2]), MOI.GreaterThan(0.0))
    vc3 = MOI.addconstraint!(m, MOI.SingleVariable(vars1[3]), MOI.GreaterThan(0.0))

    # integCons = MOI.addconstraint!(m, MOI.SingleVariable(vars1[3]), MOI.ZeroOne())
    # integCons = MOI.addconstraint!(m, MOI.SingleVariable(vars1[2]), MOI.Integer())


    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0, 5.0], vars1), 0.0)
    @test MOI.canset(m, MOI.ObjectiveFunction{typeof(saf)}())
    MOI.set!(m, MOI.ObjectiveFunction{typeof(saf)}(), saf)
    @test MOI.get(m, MOIU.AttributeFromModelCache(MOI.ObjectiveFunction{typeof(saf)}())) ≈ saf
    @test MOI.get(m, MOI.ObjectiveFunction{typeof(saf)}()) ≈ saf

    solver = CbcOptimizer()

    MOI.copy!(solver, m)
    MOI.set!(solver, MOI.ObjectiveSense(), MOI.MaxSense)
    @test MOI.get(solver, MOI.ObjectiveSense()) == MOI.MaxSense
    myStatus = MOI.optimize!(solver)

    @test MOI.get(solver, MOI.ObjectiveValue()) == 5.0
    @test MOI.get(solver, MOI.ObjectiveBound()) == 5.0
    @test MOI.get(solver, MOI.VariablePrimal(), vars1[3]) == 1

end




#
