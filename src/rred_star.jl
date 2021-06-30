"""
Return the star product as defined by Redheffer of two 2x2 BlockMatrices

See: https://en.wikipedia.org/wiki/Redheffer_star_product
"""
function rred_star(A::T, B::T) where T<:BlockMatrix
    # The inputs at the very least need to have conformable blocks
    # However, the Redheffer star is defined for blocks of all same shape
    I_AB = inv(I - A[Block(1, 2)] * B[Block(2, 1)])
    I_BA = inv(I - B[Block(2, 1)] * A[Block(1, 2)])
    return mortar(reshape([
        B[Block(1, 1)] * I_AB * A[Block(1, 1)],
        A[Block(2, 1)] + A[Block(2, 2)] * I_BA * B[Block(2, 1)] * A[Block(1, 1)],
        B[Block(1, 2)] + B[Block(1, 1)] * I_AB * A[Block(1, 2)] * B[Block(2, 2)],
        A[Block(2, 2)] * I_BA * B[Block(2, 2)]
    ], 2, 2))
end
