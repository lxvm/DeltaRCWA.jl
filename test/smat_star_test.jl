"Calculate the inverse of the S-matrix star"
function inv_smat_star(S)
    M = inv(S)
    mortar(
        (M[Block(2, 2)], M[Block(2, 1)]),
        (M[Block(1, 2)], M[Block(1, 1)]),
    )
#    mortar(
#        (inv(S[Block(2,2)] - S[Block(2,1)] * inv(S[Block(1,1)]) * S[Block(1,2)]), inv(S[Block(1,2)] - S[Block(1,1)] * inv(S[Block(2,1)]) * S[Block(2,2)])),
#        (inv(S[Block(2,1)] - S[Block(2,2)] * inv(S[Block(1,2)]) * S[Block(1,1)]), inv(S[Block(1,1)] - S[Block(1,2)] * inv(S[Block(2,2)]) * S[Block(2,1)])),
#    )
end

for n in 2:5, m in 2:5
    # setup rectangular matrices
    S = [BlockMatrix(rand(n + m, n + m), [n, m], [n, m]) for _ in 1:3]
    Lid = mortar(reshape([Matrix(a*I, n, n) for a in [0., 1., 1., 0.]], (2, 2)))
    Rid = mortar(reshape([Matrix(a*I, m, m) for a in [0., 1., 1., 0.]], (2, 2)))
    @testset "verify identity" for M in S
        @test M == smat_star(M, Rid) == smat_star(Lid, M)
    end
    @testset "verify associativity" for A in S, B in S, C in S
        # for dimension compatibility for star product of non-square blocks
        D = Array(transpose(reshape(B[end:-1:1], size(B))))
        D = BlockMatrix(D, [m, n], [m, n])
        @test smat_star(A, smat_star(D, C)) ≈ smat_star(smat_star(A, D), C)
    end
    if n == m
        @testset "verify inverse" for M in S
            M⁻¹ = inv_smat_star(M)
            @test smat_star(M, M⁻¹) ≈ smat_star(M⁻¹, M) ≈ Lid ≈ Rid
        end
    end
end
