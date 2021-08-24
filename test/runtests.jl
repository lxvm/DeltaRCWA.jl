using LinearAlgebra: I, inv
using Test: @test, @testset

using BlockArrays: mortar, Block, BlockMatrix
using LinearMaps: LinearMap

using DeltaRCWA

@testset "DeltaRCWA.jl" begin
    @testset "S-matrix  star properties" begin include("smat_star_test.jl") end
    # @testset "Redheffer star properties" begin include("rred_star_test.jl") end
end
