export smat_star, rred_star

"""
    smat_star(A::S where S <: AbstractBlockMatrix, B::T where T<: AbstractBlockMatrix)

Return the star product of two S-matrices (2x2 BlockMatrices)

See: https://en.wikipedia.org/wiki/Redheffer_star_product
"""
function smat_star(A::S where S <: AbstractBlockMatrix, B::T where T<: AbstractBlockMatrix)
    # The inputs at the very least need to have conformable blocks
    # However, the Redheffer star is defined for blocks of all same shape
    I_AB = inv(I - A[Block(2, 2)] * B[Block(1, 1)])
    I_BA = inv(I - B[Block(1, 1)] * A[Block(2, 2)])
    mortar(
        (A[Block(1, 1)] + A[Block(1, 2)] * I_BA * B[Block(1, 1)] * A[Block(2, 1)],        A[Block(1, 2)] * I_BA * B[Block(1, 2)]),
        (B[Block(2, 1)] * I_AB * A[Block(2, 1)],        B[Block(2, 2)] + B[Block(2, 1)] * I_AB * A[Block(2, 2)] * B[Block(1, 2)]),
    )
end

"""
    rred_star(A::S where S <: AbstractBlockMatrix, B::T where T<: AbstractBlockMatrix)

Return the star product as defined by Redheffer of two 2x2 BlockMatrices

See: https://en.wikipedia.org/wiki/Redheffer_star_product
"""
function rred_star(A::S where S <: AbstractBlockMatrix, B::T where T<: AbstractBlockMatrix)
    # The inputs at the very least need to have conformable blocks
    # However, the Redheffer star is defined for blocks of all same shape
    I_AB = inv(I - A[Block(1, 2)] * B[Block(2, 1)])
    I_BA = inv(I - B[Block(2, 1)] * A[Block(1, 2)])
    mortar(
        (B[Block(1, 1)] * I_AB * A[Block(1, 1)],        B[Block(1, 2)] + B[Block(1, 1)] * I_AB * A[Block(1, 2)] * B[Block(2, 2)]),
        (A[Block(2, 1)] + A[Block(2, 2)] * I_BA * B[Block(2, 1)] * A[Block(1, 1)],        A[Block(2, 2)] * I_BA * B[Block(2, 2)]),
    )
end