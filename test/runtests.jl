using Test

using deltaRCWA

@testset "deltaRCWA.jl" begin
    @testset "S-matrix  star properties" begin include("smat_star_test.jl") end
    @testset "Redheffer star properties" begin include("rred_star_test.jl") end
end
