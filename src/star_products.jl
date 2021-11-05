export ⋆ₛ, ⋆, UniformScalingExchange, J

"""
    ⋆ₛ(A::AbstractBlockMatrix, B::AbstractBlockMatrix)

Return the star product of two S-matrices (2x2 BlockMatrices)
See:
https://en.wikipedia.org/wiki/Redheffer_star_product#Connection_to_scattering_matrices
"""
function ⋆ₛ(A::AbstractBlockMatrix, B::AbstractBlockMatrix)
    # The inputs at the very least need to have conformable blocks
    # However, the Redheffer star is defined for blocks of all same shape
    invI_A₂₂B₁₁ = inv(I - A[Block(2, 2)] * B[Block(1, 1)])
    invI_B₁₁A₂₂ = inv(I - B[Block(1, 1)] * A[Block(2, 2)])
    mortar(
        (A[Block(1, 1)] + A[Block(1, 2)] * invI_B₁₁A₂₂ * B[Block(1, 1)] * A[Block(2, 1)],        A[Block(1, 2)] * invI_B₁₁A₂₂ * B[Block(1, 2)]),
        (B[Block(2, 1)] * invI_A₂₂B₁₁ * A[Block(2, 1)],        B[Block(2, 2)] + B[Block(2, 1)] * invI_A₂₂B₁₁ * A[Block(2, 2)] * B[Block(1, 2)]),
    )
end

"""
    ⋆ₛ(A::BlockMap, B::BlockMap)::BlockMap

Perform the star product on BlockMaps
"""
function ⋆ₛ(A::BlockMap, B::BlockMap)::BlockMap
    @assert (2, 2) == A.rows == B.rows "Both maps must be 2x2 block maps"
    @assert all(has_regular_blocks.([A, B])) "Both maps must have regular block alignment"
    A₁₁, A₁₂, A₂₁, A₂₂ = A.maps
    B₁₁, B₁₂, B₂₁, B₂₂ = B.maps
    invI_A₂₂B₁₁ = LinearMap(x -> gmres(I - A₂₂ * B₁₁, x), size(A₂₂, 1), size(B₁₁, 2))
    invI_B₁₁A₂₂ = LinearMap(x -> gmres(I - B₁₁ * A₂₂, x), size(B₁₁, 1), size(A₂₂, 2))
    [
        A₁₁ + A₁₂ * invI_B₁₁A₂₂ * B₁₁ * A₂₁     A₁₂ * invI_B₁₁A₂₂ * B₁₂;
        B₂₁ * invI_A₂₂B₁₁ * A₂₁                 B₂₂ + B₂₁ * invI_A₂₂B₁₁ * A₂₂ * B₁₂
    ]
end

"""
    has_regular_blocks(M::BlockMap)::Bool

Returns true if the blocks in the BlockMap are partitioned by a row partition and a
column partition (i.e. all blocks in each column/row have the same width/height).
"""
function has_regular_blocks(M::BlockMap)::Bool
    Ncol = M.rows[1]
    Nrow = length(M.rows)
    all(Ncol == M.rows[i] for i in 1:Nrow) &&
    all(M.colranges[i] == M.colranges[i+j*Ncol] for i in 1:Ncol, j in 0:(Nrow-1))
end

"""
    ⋆ₛ(A::LinearMap, B::LinearMap)::LinearMap

Take the scattering star product of two matrix-free scattering matrices.
Requires the inputs to be square matrices of the same size and assumes that the
matrix is divided into equal-size quadrants to obtain the blocks
"""
function ⋆ₛ(A::LinearMap, B::LinearMap)::LinearMap
    d = Int(size(A, 1) // 2)
    @assert size(A) == size(B) == reverse(size(A)) "maps must be endomorphisms"
    A₁₁ = LinearMap(x -> (A * vcat(x, zero(x)))[1:d], d)
    A₁₂ = LinearMap(x -> (A * vcat(zero(x), x))[1:d], d)
    A₂₂ = LinearMap(x -> (A * vcat(zero(x), x))[(1+d):2d], d)
    B₁₁ = LinearMap(x -> (B * vcat(x, zero(x)))[1:d], d)
    B₂₁ = LinearMap(x -> (B * vcat(x, zero(x)))[(1+d):2d], d)
    B₂₂ = LinearMap(x -> (B * vcat(zero(x), x))[(1+d):2d], d)
    LinearMap(2d) do x
        x₁ = @view x[1:d]
        x₂ = @view x[(d+1):2d]
        # This is slightly optimized to avoid repeating two linear mappings
        y = A * vcat(x₁, fill(zero(eltype(x)), d))
        y₁ = @view y[1:d]
        y₂ =  gmres(I - A₂₂ * B₁₁, y[(d+1):2d])
        z = B * vcat(fill(zero(eltype(x)), d), x₂)
        z₁ = gmres(I - B₁₁ * A₂₂, z[1:d])
        z₂ = @view z[(d+1):2d]
        vcat(
            y₁ + A₁₂ * B₁₁ * y₂ + A₁₂ * z₁,
            B₂₁ * y₂ + z₂ + B₂₁ * A₂₂ * z₁
        )
    end
end

# """
#     exchangehalves(A::AbstractArray{<:Any,N}, dim::Integer) where N

# Exchange the first half of A along first dimension with the second half.
# Will fail if size(A, dim) is odd.
# """
# function exchangehalves(A::AbstractArray{<:Any, N}, dim::Integer) where N
#     if N == 0
#         return A
#     else
#         @assert 1 <= dim <= N
#         return _exchangehalves(A, dim)
#     end
# end
# function _exchangehalves(A::AbstractArray, dim)
#     d = Int(size(A, dim) // 2)
#     selectdim(A, dim, vcat((d+1):2d, 1:d))
# end
# # function _exchangehalves(A::BlockArray, dim)
# #     @assert length(blockaxes(A, dim)) == 2
# #     selectdim(A, dim, blockslice)
# # end

# """
#     ExchangeOperator

# Type equivalent to the permutation matrix

#     [
#         0 1;
#         1 0
#     ]

# with flexible size set so that each block is the same square size.
# """
# struct ExchangeOperator end

# """
#     J

# Object of type ExchangeOperator
# """
# const J = ExchangeOperator()

# """
#     ExchangeOperator()(d::Int)

# Returns a LinearMap version of the exchangehalves operator. `d` is the axis length
# and must be even.
# """
# (::ExchangeOperator, d::Int) = LinearMap{Bool}(d; issymmetric=true, ishermitian=true) do A
#     @assert size(A, 1) == d
#     exchangehalves(A, 1)
# end

# Base.:*(::ExchangeOperator, A::AbstractArray) = exchangehalves(A, 1)
# Base.:*(A::AbstractArray{<:Any,N}, ::ExchangeOperator) where N = exchangehalves(A, N)

# Base.:*(::ExchangeOperator, A::LinearMap) = J(size(A, 1)) * A
# Base.:*(A::LinearMap, ::ExchangeOperator) = A * J(size(A, 2))

# """
#     ⋆(A, B)

# Return the star product of two matrices
# See:
# https://en.wikipedia.org/wiki/Redheffer_star_product
# """
# ⋆(A, B) = J * ((J * A) ⋆ₛ (J * B))