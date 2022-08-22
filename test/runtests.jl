#
using Test, SafeTestsets

@testset "NodalPolynomialSpaces.jl" begin

@time @safetestset "Interface" begin include("interface.jl") end
@time @safetestset "Examples" begin include("examples.jl") end

end
