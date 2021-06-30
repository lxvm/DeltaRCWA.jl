for n in 4:10
    # setup square matrices
    S = [BlockMatrix(rand(2n, 2n), [n, n], [n, n]) for _ in 1:3]
    id = BlockMatrix(Matrix(1.0I, 2n, 2n), [n, n], [n, n])
    @testset "verify identity" for M in S
        @test M == rred_star(M, id) == rred_star(id, M)
    end
    @testset "verify associativity" for A in S, B in S, C in S
        @test (rred_star(A, rred_star(B, C)) ≈ rred_star(rred_star(A, B), C))
    end
    @testset "verify inverse" for M in S
        M⁻¹ = inv(M)
        @test rred_star(M, M⁻¹) ≈ rred_star(M⁻¹, M) ≈ id
    end
end

