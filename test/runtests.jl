using Cbc, Compat.Test

@testset "MPB" begin
    include("MPB_wrapper.jl")
end

@testset "MOI" begin
    include("MOI_wrapper.jl")
end
