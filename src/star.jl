using LinearAlgebra

using BlockArrays

function star(A::T, B::T) where T<:BlockMatrix
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

# Example matrix A of size (2N, 2N)
# mortar(reshape([rand(Float64, (N, N)) for i in 1:4], 2, 2))

# Example identity
# BlockMatrix(Matrix(Diagonal(fill(1., 2N))), [N, N], [N, N])

# Tests (to be completed)
# identity
# associativity
# inverse

# Idea for numeric types: we have type inclusions Integer<:Real<:Number
# but Integers for a ring and Rational, Real, Complex numbers are fields
# so maybe Real<:Union{Number, Field} and so on to reflect the mathematical
# structure of these types (but of course a machine is only finite)
# For instance, a field guarantees inverses in vector spaces (but not Ints)
# I have no idea whether this is convenient
# Related:
# https://github.com/JuliaLang/julia/issues/5
# https://discourse.julialang.org/t/struct-as-subtype-of-multiple-types/52357/2
# https://discourse.julialang.org/t/parametric-abstract-types/25651/2
