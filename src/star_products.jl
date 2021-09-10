export smat_star, ⋆


⋆ = function smat_star end

"""
    smat_star(A::AbstractBlockMatrix, B::AbstractBlockMatrix)

Return the star product of two S-matrices (2x2 BlockMatrices)
See:
https://en.wikipedia.org/wiki/Redheffer_star_product#Connection_to_scattering_matrices
"""
function smat_star(A::AbstractBlockMatrix, B::AbstractBlockMatrix)
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
    smat_star(A::BlockMap, B::BlockMap)::BlockMap

Perform the star product on BlockMaps
"""
function smat_star(A::BlockMap, B::BlockMap)::BlockMap
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
    smat_star(A::LinearMap, B::LinearMap)::LinearMap

Take the scattering star product of two matrix-free scattering matrices.
Requires the inputs to be square matrices of the same size and assumes that the
matrix is divided into equal-size quadrants to obtain the blocks
"""
function smat_star(A::LinearMap, B::LinearMap)::LinearMap
    d = Int(size(A, 1) // 2)
    @assert size(A) == size(B) == reverse(size(A)) "maps must be endomorphisms"
    A₁₁ = LinearMap(x -> (A * vcat(x, zero(x)))[1:d], d)
    A₁₂ = LinearMap(x -> (A * vcat(zero(x), x))[1:d], d)
    A₂₂ = LinearMap(x -> (A * vcat(zero(x), x))[(1+d):2d], d)
    B₁₁ = LinearMap(x -> (B * vcat(x, zero(x)))[1:d], d)
    B₂₁ = LinearMap(x -> (B * vcat(x, zero(x)))[(1+d):2d], d)
    B₂₂ = LinearMap(x -> (B * vcat(zero(x), x))[(1+d):2d], d)
    invI_A₂₂B₁₁ = LinearMap(x -> gmres(I - A₂₂ * B₁₁, x), size(A₂₂, 1), size(B₁₁, 2))
    invI_B₁₁A₂₂ = LinearMap(x -> gmres(I - B₁₁ * A₂₂, x), size(B₁₁, 1), size(A₂₂, 2))
    LinearMap(2d) do x
        x₁ = @view x[1:d]
        x₂ = @view x[(d+1):2d]
        # This is slightly optimized to avoid repeating two linear mappings
        y = A * vcat(x₁, fill(zero(eltype(x)), d))
        y₁ = @view y[1:d]
        y₂ = @view y[(d+1):2d]
        z = B * vcat(fill(zero(eltype(x)), d), x₂)
        z₁ = @view z[1:d]
        z₂ = @view z[(d+1):2d]
        vcat(
            y₁ + A₁₂ * invI_B₁₁A₂₂ * B₁₁ * y₂ + A₁₂ * invI_B₁₁A₂₂ * z₁,
            B₂₁ * invI_A₂₂B₁₁ * y₂ + z₂ + B₂₁ * invI_A₂₂B₁₁ * A₂₂ * z₁
        )
    end
end