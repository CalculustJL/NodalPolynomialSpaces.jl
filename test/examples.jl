#
@testset "1D tests" begin
    @time @safetestset "Heat" begin include("../examples/1D_Poisson.jl") end
end

@testset "2D tests" begin
    @time @safetestset "Heat trans" begin include("../examples/2D_Poisson.jl") end
end
#
