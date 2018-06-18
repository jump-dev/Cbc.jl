using Base.Test, MathOptInterface, MathOptInterface.Test, MathOptInterface.Utilities

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test
const MOIU = MathOptInterface.Utilities
const MOIB = MathOptInterface.Bridges

@MOIU.model ModelForCachingOptimizer (ZeroOne, Integer) (EqualTo, GreaterThan, LessThan, Interval) () () (SingleVariable,) (ScalarAffineFunction,) () ()

const optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), CbcOptimizer())
const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals = false, infeas_certificates = false)

@testset "Continuous linear problems" begin
    # AlmostSuccess for linear9 with SCS 2
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config)
end





#
# using Base.Test, MathOptInterface, MathOptInterface.Test, MathOptInterface.Utilities
#
# const MOI  = MathOptInterface
# const MOIT = MathOptInterface.Test
# const MOIU = MathOptInterface.Utilities
#
# @MOIU.model ModelForCachingOptimizer (ZeroOne, Integer) (EqualTo, GreaterThan, LessThan, Interval) () () (SingleVariable,) (ScalarAffineFunction,) () ()
#
# const optimizer = MOIU.CachingOptimizer(ModelForCachingOptimizer{Float64}(), CbcOptimizer())
# const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals = false, infeas_certificates = false)
#
# @testset "Contlinear tests" begin
#     # @testset "linear1test" begin MOIT.linear1test(optimizer, config) end
#     # @testset "linear2test" begin MOIT.linear2test(optimizer, config) end
#     # @testset "linear3test" begin MOIT.linear3test(optimizer, config) end
#     # @testset "linear4test" begin MOIT.linear4test(optimizer, config) end
#     # @testset "linear5test" begin MOIT.linear5test(optimizer, config) end
#     # @testset "linear6test" begin MOIT.linear6test(optimizer, config) end
#     # @testset "linear8atest" begin MOIT.linear8atest(optimizer, config) end
#     # @testset "linear8btest" begin MOIT.linear8btest(optimizer, config) end
#     # @testset "linear8ctest" begin MOIT.linear8ctest(optimizer, config) end
#     # @testset "linear9test" begin MOIT.linear9test(optimizer, config) end
#     # @testset "linear10test" begin MOIT.linear10test(optimizer, config) end
#     @testset "linear11test" begin MOIT.linear11test(optimizer, config) end
#     # @testset "linear12test" begin MOIT.linear12test(optimizer, config) end
#     @testset "linear13test" begin MOIT.linear13test(optimizer, config) end
#     @testset "linear14test" begin MOIT.linear14test(optimizer, config) end
# end


#
