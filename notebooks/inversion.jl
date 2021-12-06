### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 16ff4dea-b7e6-4165-adaf-b70754af7800
using LinearAlgebra

# ╔═╡ 9643b8cc-2fc9-462f-9614-c92c71ee9c14
using BlockArrays

# ╔═╡ 7a4eed4d-e1e1-4bba-aa96-ee2de4812352
using KrylovKit

# ╔═╡ 7c81445c-a46b-4d38-82a9-54f83d9ef539
using IterativeSolvers

# ╔═╡ 38cfc033-dcbd-4e92-832e-5066e095087c
using LinearMaps

# ╔═╡ 71d5a8a5-d4d2-4261-b0fd-fa6a962a9315
using Plots

# ╔═╡ 47f862e0-31a1-4448-93a8-8fa0d990f5ba
using BenchmarkTools

# ╔═╡ 4943174a-004c-11ec-3f4b-398aac2604a4
md"
# Inversion in Redheffer star products

The goal of this notebook is to verify tricks for inverting blocks of matrices
using matrix-free iterative solvers for use in computing Redheffer star products.
"

# ╔═╡ ca594862-6328-4d68-a80f-dbd6a76742f1
md"
## The problem

Given a scattering matrix 

$S = 
\begin{pmatrix}
S_{11} & S_{12}
\\
S_{21} & S_{22}
\end{pmatrix}$

and a [star product in the scattering matrix convention](
https://en.wikipedia.org/wiki/Redheffer_star_product#Connection_to_scattering_matrices)

$A \star_S B =
\begin{pmatrix}
    A_{11} + A_{12} (I - B_{11} A_{22})^{-1} B_{11} A_{21} &
    A_{12} (I - B_{11} A_{22})^{-1} B_{12}
    \\
    B_{21} (I - A_{22} B_{11})^{-1} A_{21} &
    B_{22} + B_{21} (I - A_{22} B_{11})^{-1} A_{22} B_{12}
\end{pmatrix}$

I am working with scattering matrices of the form $S = A^{-1} B$, so in practice
to compute star products I will need to be able to apply the linear transformation
in each block of $A^{-1}B$.
This form arises from GSTCs.
For propagation matrices or for scattering at a dielectric interface, $S$ is known.
"

# ╔═╡ e2772ee2-2b2a-4fba-8a2a-163a14e2e00f
md"
## Closed form block inversions

[Redheffer published a formula](https://doi.org/10.1002%2Fsapm1960391269)
for the matrix inverse for a 2x2 block matrix

$A^{-1} = 
\begin{pmatrix}
    (A_{11} - A_{12} A_{22}^{-1} A_{21})^{-1}
	&
    (A_{21} - A_{22} A_{12}^{-1} A_{11})^{-1}
    \\
    (A_{12} - A_{11} A_{21}^{-1} A_{22})^{-1}
	&
    (A_{22} - A_{21} A_{11}^{-1} A_{12})^{-1}
\end{pmatrix}$

In proving that this matrix is indeed the inverse, one may find it useful to know the
[Woodbury Matrix Identity](https://en.wikipedia.org/wiki/Woodbury_matrix_identity).

This formula is useful to test against, however it would be impractical to implement
with iterative solvers due to the nested inversions making the complexity worse.
Since our goal is to have a matrix-free formulation, we will only test against this
formula for dense matrices.

For a generic inversion algorithm, such as
[GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method),
the [time complexity is quadratic in the number of iterations](
https://scicomp.stackexchange.com/questions/30681/operation-count-for-gmres),
which is what makes nested inversions more costly.
(For symmetric matrices (not my case), consider Conjugate-Gradient algorithms.)
"

# ╔═╡ 1d8d4f1f-a80f-4708-add5-bc7940f068ba
md"
## Iterative block inversions

Steven observed that to calculate the matrix-vector product with a block of $A$
one can get away with using iterative solvers to invert $A$, which should be fast.

The process is as follows: suppose we wish to calculate

$y_1 = (A^{-1})_{11} x_1$
$y_2 = (A^{-1})_{21} x_1$

then instead solve the inversion problem

$\begin{pmatrix}
y_1
\\
y_2
\end{pmatrix}
= A^{-1}
\begin{pmatrix}
x_1
\\
0
\end{pmatrix}$

Similarly solve

$z_1 = (A^{-1})_{12} x_2$
$z_2 = (A^{-1})_{22} x_2$

with

$\begin{pmatrix}
z_1
\\
z_2
\end{pmatrix}
= A^{-1}
\begin{pmatrix}
0
\\
x_2
\end{pmatrix}$
"

# ╔═╡ 06b86080-af86-426d-8777-b3326c19b6aa
md"
## Comparisons

We will now compare the results of the methods discussed above for accuracy.

Note that for many iterative methods to work, the eigenvalues should not be clustered
at the origin because the convergence will be slow (the power method won't work /  
Krylov subspace will be influenced by too many eigenvalues of similar norm).
[Explanation in here](https://courses.maths.ox.ac.uk/node/view_material/53060).
"

# ╔═╡ dac13df7-27e3-4367-8db8-bb8cb7439364
begin
	# setup arbitrary linear transformation
	n = 50
	### test with a random matrix and translate the characteristic polynomial
	A = rand(2n, 2n) + 10I(2n)
	### symmetrize matrix
	# A = transpose(A) * A
	# extract blocks
	A₁₁ = A[1:n, 1:n]
	A₁₂ = A[1:n, (n+1):end]
	A₂₁ = A[(n+1):end, 1:n]
	A₂₂ = A[(n+1):end, (n+1):end]
end;

# ╔═╡ 8a53947d-87f0-4163-9c50-32d23a5d144b
plot(eigvals(A), seriestype=:scatter)

# ╔═╡ 238d9089-9dda-4750-874b-502e11c35138
begin
	# compute inverted blocks using dense methods
	invA = inv(A)
	invA₁₁ = invA[1:n, 1:n]
	invA₁₂ = invA[1:n, (n+1):end]
	invA₂₁ = invA[(n+1):end, 1:n]
	invA₂₂ = invA[(n+1):end, (n+1):end]
end;

# ╔═╡ 0d19e309-c54f-4e66-b172-142aa2ffb29f
begin
	# compute inverted blocks using Redheffer's formula
	RinvA₁₁ = inv(A₁₁ - A₁₂ * inv(A₂₂) * A₂₁)
	RinvA₁₂ = inv(A₂₁ - A₂₂ * inv(A₁₂) * A₁₁)
	RinvA₂₁ = inv(A₁₂ - A₁₁ * inv(A₂₁) * A₂₂)
	RinvA₂₂ = inv(A₂₂ - A₂₁ * inv(A₁₁) * A₁₂)
end;

# ╔═╡ 61b00eb3-fca2-46a2-8ddb-080c819a4670
# verify Redheffer and dense methods match
all([
	invA₁₁ ≈ RinvA₁₁,
	invA₁₂ ≈ RinvA₁₂,
	invA₂₁ ≈ RinvA₂₁,
	invA₂₂ ≈ RinvA₂₂,
])

# ╔═╡ eb9cf21e-7bda-49be-847a-7737b84638a4
function onehotvector(i::T, n::T, ::S=true) where {T <: Integer, S}
	out = zeros(S, n)
	out[i] = one(S)
	out
end

# ╔═╡ 8163e2bb-b072-4284-8da9-29af6f212ba2
# verify linsolve returns columns of inverse matrix
collect(invA[:, i] ≈ linsolve(A, onehotvector(i, 2n, 1.0))[1] for i in 1:2n)

# ╔═╡ f5a0443d-fe89-41cf-913b-6e278f34d48e
# verify gmres returns columns of inverse matrix
collect(invA[:, i] ≈ gmres(A, onehotvector(i, 2n, 1.0)) for i in 1:2n)

# ╔═╡ 0eec9c1f-2546-410b-850b-c3de8ded3e30
# verify inversion against some random vector
s = randn(2n)

# ╔═╡ 088c4638-6326-4188-9827-1c62135d9c6d
A * inv(A) * s ≈ s

# ╔═╡ 39f83ea9-92ea-4b41-9d23-17ef8b74988a
x, info = linsolve(A, s)

# ╔═╡ 1c1dce9c-db71-418e-836a-3b01650c5a66
A * x ≈ s

# ╔═╡ 827df9d6-93fa-4d30-82c4-d78314b50013
info

# ╔═╡ 50749a6e-fcb3-4427-b3ee-96fe2549898a
A * gmres(A, s) ≈ s

# ╔═╡ 1648540a-6812-47c1-ab86-2d67c6fb4838
md"
## Benchmarking

We will now do a brief benchmark of the methods compared above
"

# ╔═╡ 69006ecc-d67c-4654-ab58-123e39f8d300
@benchmark inv(A)

# ╔═╡ 10919814-a0b6-423f-ab60-7baa8218f9f5
@benchmark linsolve(A, s)

# ╔═╡ 0e916044-6944-4e75-a944-a6669851ed21
@benchmark gmres(A, s)

# ╔═╡ c8de38f2-5dc4-4abd-a735-11df723e6b85
L = LinearMap(x -> A₁₁ * x + A₁₂ * gmres(A₂₂, A₂₁ * x), n)

# ╔═╡ 487efe4c-b847-4e62-a0b6-78c759b90bf4
### inverting a single block of the matrix (5-10 times slower)
@benchmark gmres(L, s[1:n])

# ╔═╡ f5a30faf-6cc6-4708-a48d-b3839f437d00
md"
## Implementing a star product

We can implement a star product with these block operations.
To do so we can use the [`LinearMaps.BlockMap`](
https://jutho.github.io/LinearMaps.jl/stable/types/#BlockMap-and-BlockDiagonalMap)
type.
"

# ╔═╡ c9f91c42-238c-4f03-8d7d-a837002357d8
M = [
	L 3L;
	2L 4L
]

# ╔═╡ 5f84cc39-1e39-4c24-95c4-50899defd70d
"""
	has_regular_blocks(M::LinearMaps.BlockMap)::Bool

Returns true if the blocks in the BlockMap are partitioned by a row partition and a
column partition (i.e. all blocks in each column/row have the same width/height).
"""
function has_regular_blocks(M::LinearMaps.BlockMap)::Bool
	Ncol = M.rows[1]
	Nrow = length(M.rows)
	all(Ncol == M.rows[i] for i in 1:Nrow) && all(M.colranges[i] == M.colranges[i+j*Ncol] for i in 1:Ncol, j in 0:(Nrow-1))
end

# ╔═╡ 56c8e6e0-d286-4210-a67e-3ddef8228dd1
begin
	nrow, ncol = 1, 3
	has_regular_blocks(hvcat(Tuple(ncol for _ in 1:nrow), fill(L, nrow*ncol)...))
end

# ╔═╡ a2500f48-2406-47f1-adec-bd6e50f67955
function ⋆(A::LinearMaps.BlockMap, B::LinearMaps.BlockMap)::LinearMaps.BlockMap
	@assert (2, 2) == A.rows == B.rows "Both maps must be 2x2 block maps"
	@assert all(has_regular_blocks.([A, B])) "Both maps must have regular block alignment"
	A₁₁, A₁₂, A₂₁, A₂₂ = A.maps
	B₁₁, B₁₂, B₂₁, B₂₂ = B.maps
	invI_A₂₂B₁₁ = LinearMap(x -> gmres(I - A₂₂ * B₁₁, x), size(A₂₂, 1), size(B₁₁, 2))
    invI_B₁₁A₂₂ = LinearMap(x -> gmres(I - B₁₁ * A₂₂, x), size(B₁₁, 1), size(A₂₂, 2))
	[
		A₁₁ + A₁₂ * invI_B₁₁A₂₂ * B₁₁ * A₂₁ 	A₁₂ * invI_B₁₁A₂₂ * B₁₂;
	 	B₂₁ * invI_A₂₂B₁₁ * A₂₁ 				B₂₂ + B₂₁ * invI_A₂₂B₁₁ * A₂₂ * B₁₂
	]
end

# ╔═╡ 1a23e5ce-c6e3-4658-83ad-9060ebe44230
id = [
	0I(n)	1I(n);
	1I(n)	LinearMap(0I(n))
]

# ╔═╡ feba7f54-b0de-4761-919d-1260dcc14b8f
md"
## Layers of nesting

We should confirm that the star products return the same result within a tolerance
when there are several (>2) matrices to combine.
To do this comparison we will define additional methods.
"

# ╔═╡ b952bfe3-2bd8-44be-a81a-b30a4f48e14a
function ⋆(A::AbstractBlockMatrix, B::AbstractBlockMatrix)
    # The inputs at the very least need to have conformable blocks
    # However, the Redheffer star is defined for blocks of all same shape
    invI_A₂₂B₁₁ = inv(I - A[Block(2, 2)] * B[Block(1, 1)])
    invI_B₁₁A₂₂ = inv(I - B[Block(1, 1)] * A[Block(2, 2)])
    mortar(
        (A[Block(1, 1)] + A[Block(1, 2)] * invI_B₁₁A₂₂ * B[Block(1, 1)] * A[Block(2, 1)],        A[Block(1, 2)] * invI_B₁₁A₂₂ * B[Block(1, 2)]),
        (B[Block(2, 1)] * invI_A₂₂B₁₁ * A[Block(2, 1)],        B[Block(2, 2)] + B[Block(2, 1)] * invI_A₂₂B₁₁ * A[Block(2, 2)] * B[Block(1, 2)]),
    )
end

# ╔═╡ d8748cf3-94df-4875-9716-44769c83f3b2
function ⋆(A::LinearMap, B::LinearMap)::LinearMap
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

# ╔═╡ 10f7b132-674e-4655-91f4-1b3953eefcdc
inv(I - A₂₂ * A₁₁) ≈ Matrix(LinearMap(size(A₁₁, 1), size(A₂₂, 2)) do x; gmres(I - A₂₂ * A₁₁, x); end)

# ╔═╡ 99420d2d-2eaf-4dd3-9dc7-ea08fa955be5
md"
check that block indexing is correct
"

# ╔═╡ 0d13c2e4-59c7-43fb-a918-0715f232a406
function ⋆(A, B)
    d = Int(size(A, 1) // 2)
    @assert size(A) == size(B) == reverse(size(A)) "maps must be endomorphisms"
    A₁₁ = LinearMap(x -> (A * vcat(x, zero(x)))[1:d], d)
    A₁₂ = LinearMap(x -> (A * vcat(zero(x), x))[1:d], d)
    A₂₂ = LinearMap(x -> (A * vcat(zero(x), x))[(1+d):2d], d)
    B₁₁ = LinearMap(x -> (B * vcat(x, zero(x)))[1:d], d)
    B₂₁ = LinearMap(x -> (B * vcat(x, zero(x)))[(1+d):2d], d)
    B₂₂ = LinearMap(x -> (B * vcat(zero(x), x))[(1+d):2d], d)
    invI_A₂₂B₁₁ = LinearMap(inv(Matrix(I - A₂₂ * B₁₁)))
    # invI_A₂₂B₁₁ = LinearMap(x -> gmres(I - A₂₂ * B₁₁, x), size(A₂₂, 1), size(B₁₁, 2))
    invI_B₁₁A₂₂ = LinearMap(inv(Matrix(I - B₁₁ * A₂₂)))
    # invI_B₁₁A₂₂ = LinearMap(x -> gmres(I - B₁₁ * A₂₂, x), size(B₁₁, 1), size(A₂₂, 2))
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

# ╔═╡ 774d6cbb-07ce-4a9c-bf82-b7fd475dcca6
# verify the star identity
M * ones(2n) ≈ (M ⋆ id) * ones(2n) ≈ (id ⋆ M) * ones(2n)

# ╔═╡ 4d292003-38af-4060-bd5c-6a6a1e2244f1
begin
	Nlayers = 2
	blocks = [Tuple(100I - rand(n, n) for _ in 1:4) for _ in 1:Nlayers]
    bmat = [mortar(e[[1, 3]], e[[2, 4]]) for e in blocks]
	bmap = [[LinearMap(e[1]) e[3]; e[2] e[4]] for e in blocks]
	lmap = [LinearMap([e[1] e[3]; e[2] e[4]]) for e in blocks]
	map = [[e[1] e[3]; e[2] e[4]] for e in blocks]
end;

# ╔═╡ 368ce9a1-4ec8-4f9c-a148-68325a4b7194
blocks[1][1] ≈ Matrix(LinearMap(n,n) do x; (map[1] * vcat(x, zero(x)))[1:n]; end)

# ╔═╡ 18ad0ba7-0d61-49c2-b9eb-9ef2a739e40d
blocks[1][2] ≈ Matrix(LinearMap(n,n) do x; (map[1] * vcat(x, zero(x)))[(1+n):2n]; end)

# ╔═╡ 903130b1-980b-4259-868a-66f207d48e41
blocks[1][3] ≈ Matrix(LinearMap(n,n) do x; (map[1] * vcat(zero(x), x))[1:n]; end)

# ╔═╡ 47c10091-3b39-4f6a-bdde-0b1bbcf10e30
blocks[1][4] ≈ Matrix(LinearMap(n,n) do x; (map[1] * vcat(zero(x), x))[(1+n):2n]; end)

# ╔═╡ cc9b6619-1f7d-462e-962d-5027eb8b5c3c
# for el in (:bmat, :bmap, :lmap)
# 	@eval $(Symbol(:S, el)) = $(el)[1]
# 	for j in 2:Nlayers
# 		@eval $(Symbol(:S, el)) = $(Symbol(:S, el)) ⋆ $(el)[$j]
# 	end
# end

# ╔═╡ 823dd2f5-d3bd-4069-bbe1-ad838ddca4a8
begin
	Sbmat = bmat[1]
	for j in 2:Nlayers
		Sbmat = Sbmat ⋆ bmat[j]
	end
end

# ╔═╡ bd0905fe-00bf-4ad0-b0d1-5d5d9aba04cf
begin
	Sbmap = bmap[1]
	for j in 2:Nlayers
		Sbmap = Sbmap ⋆ bmap[j]
	end
end

# ╔═╡ 2a06b4e8-b3e8-436f-848b-bde664175505
begin
	Slmap = lmap[1]
	for j in 2:Nlayers
		Slmap = Slmap ⋆ lmap[j]
	end
end

# ╔═╡ 4b9000b8-2210-48ae-afef-08055fa1851d
begin
	Smap = map[1]
	for j in 2:Nlayers
		Smap = Smap ⋆ map[j]
	end
end

# ╔═╡ 374a89d2-708f-4d6b-9daa-49363a68f218
v = rand(2n)

# ╔═╡ 5a9d5d42-ddda-4b2e-aeae-020a39ae7c35
vbmat = Sbmat * v

# ╔═╡ 14a0a868-c6ca-40a6-95fd-4e2052936ca3
vbmap = Sbmap * v

# ╔═╡ a5aadcc0-8b47-4dde-8373-800d4583e369
vlmap = Slmap * v

# ╔═╡ 3d058953-c4a8-453b-b0ca-7ad8fa930c48
vmap = Smap * v

# ╔═╡ 0ad36213-5639-443e-81d3-d8904f507639
Smap

# ╔═╡ e501ad52-3ed5-4151-8dc1-7a362158e8c3
vbmat ≈ vbmap

# ╔═╡ 3b2038e7-87a7-4bd7-b823-f230058c95b9
vbmat ≈ vlmap

# ╔═╡ eb74c6e9-1db0-4dbc-8222-47b5fc8744db
vbmat ≈ vmap

# ╔═╡ 96f6d024-5c55-4b3b-adbf-4cc8eb996adb
vlmap ≈ vmap

# ╔═╡ e0c2c9a3-84be-4857-8e2a-022e29a1ae77
vbmap ≈ vmap

# ╔═╡ c4b7d46c-411c-462b-b326-d2816c61c8c9
vbmap ≈ vlmap

# ╔═╡ 4b3a570c-ba16-4bda-b8d8-81ade4b5e3e7
norm(vmap - vlmap) / (norm(vmap) * norm(vlmap))

# ╔═╡ cc4fcea9-eefe-4f5c-95bc-a1eac6873767
begin
	plot(vbmat, label="blockarray")
	plot!(vbmap, label="block map")
	plot!(vlmap, label="linear map")
	plot!(vmap, label="array")
end

# ╔═╡ d25871e5-5e97-4b56-9c1c-0eac45d778d2
md"
## More facts

Similarly, we can compute the star inverse

$S^{-\star} \star S = S \star S^{-\star} = I_\star$

where

$I_\star = 
\begin{pmatrix}
0 & I
\\
I & 0
\end{pmatrix}$

and

$S^{-\star} = 
\begin{pmatrix}
    (S_{22} - S_{21} S_{11}^{-1} S_{12})^{-1}
	&
    (S_{12} - S_{11} S_{21}^{-1} S_{22})^{-1}
    \\
    (S_{21} - S_{22} S_{12}^{-1} S_{11})^{-1}
	&
    (S_{11} - S_{12} S_{22}^{-1} S_{21})^{-1}
\end{pmatrix}$
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearMaps = "7a12625a-238d-50fd-b39a-03d52299707e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
BenchmarkTools = "~1.1.1"
BlockArrays = "~0.15.3"
IterativeSolvers = "~0.9.1"
KrylovKit = "~0.5.3"
LinearMaps = "~3.4.0"
Plots = "~1.20.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "623a32b87ef0b85d26320a8cc7e57ded707aef64"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.7.5"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "42ac5e523869a84eac9669eaceed9e4aa0e1587b"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.1.4"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "a37151e369c618aebaff8b95b9db2f603246e160"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "0.15.3"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "693210145367e7685d8604aee33d9bfb85db8b31"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.11.9"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "182da592436e287758ded5be6e32c406de3a2e47"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "74ef6288d071f58033d54fd6708d4bc23a8b8972"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "0328ad9966ae29ccefb4e1b9bfd8c8867e4360df"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearMaps]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "a03a9d24d4597725d9e826cee1ef5164e6bdd661"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.4.1"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "8365fa7758e2e8e4443ce866d6106d8ecbb4474e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.20.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "1f27772b89958deed68d2709e5f08a5e5f59a5af"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.7"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "0f2aa8e32d511f758a2ce49208181f7733a0936a"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.1.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2bb0cb32026a66037360606510fca5984ccc6b75"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.13"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─4943174a-004c-11ec-3f4b-398aac2604a4
# ╠═16ff4dea-b7e6-4165-adaf-b70754af7800
# ╠═9643b8cc-2fc9-462f-9614-c92c71ee9c14
# ╠═7a4eed4d-e1e1-4bba-aa96-ee2de4812352
# ╠═7c81445c-a46b-4d38-82a9-54f83d9ef539
# ╠═38cfc033-dcbd-4e92-832e-5066e095087c
# ╠═71d5a8a5-d4d2-4261-b0fd-fa6a962a9315
# ╟─ca594862-6328-4d68-a80f-dbd6a76742f1
# ╟─e2772ee2-2b2a-4fba-8a2a-163a14e2e00f
# ╟─1d8d4f1f-a80f-4708-add5-bc7940f068ba
# ╟─06b86080-af86-426d-8777-b3326c19b6aa
# ╠═dac13df7-27e3-4367-8db8-bb8cb7439364
# ╠═8a53947d-87f0-4163-9c50-32d23a5d144b
# ╠═238d9089-9dda-4750-874b-502e11c35138
# ╠═0d19e309-c54f-4e66-b172-142aa2ffb29f
# ╠═61b00eb3-fca2-46a2-8ddb-080c819a4670
# ╠═eb9cf21e-7bda-49be-847a-7737b84638a4
# ╠═8163e2bb-b072-4284-8da9-29af6f212ba2
# ╠═f5a0443d-fe89-41cf-913b-6e278f34d48e
# ╠═0eec9c1f-2546-410b-850b-c3de8ded3e30
# ╠═088c4638-6326-4188-9827-1c62135d9c6d
# ╠═39f83ea9-92ea-4b41-9d23-17ef8b74988a
# ╠═1c1dce9c-db71-418e-836a-3b01650c5a66
# ╠═827df9d6-93fa-4d30-82c4-d78314b50013
# ╠═50749a6e-fcb3-4427-b3ee-96fe2549898a
# ╟─1648540a-6812-47c1-ab86-2d67c6fb4838
# ╠═47f862e0-31a1-4448-93a8-8fa0d990f5ba
# ╠═69006ecc-d67c-4654-ab58-123e39f8d300
# ╠═10919814-a0b6-423f-ab60-7baa8218f9f5
# ╠═0e916044-6944-4e75-a944-a6669851ed21
# ╠═c8de38f2-5dc4-4abd-a735-11df723e6b85
# ╠═487efe4c-b847-4e62-a0b6-78c759b90bf4
# ╟─f5a30faf-6cc6-4708-a48d-b3839f437d00
# ╠═c9f91c42-238c-4f03-8d7d-a837002357d8
# ╠═5f84cc39-1e39-4c24-95c4-50899defd70d
# ╠═56c8e6e0-d286-4210-a67e-3ddef8228dd1
# ╠═a2500f48-2406-47f1-adec-bd6e50f67955
# ╠═1a23e5ce-c6e3-4658-83ad-9060ebe44230
# ╠═774d6cbb-07ce-4a9c-bf82-b7fd475dcca6
# ╟─feba7f54-b0de-4761-919d-1260dcc14b8f
# ╠═b952bfe3-2bd8-44be-a81a-b30a4f48e14a
# ╠═d8748cf3-94df-4875-9716-44769c83f3b2
# ╠═10f7b132-674e-4655-91f4-1b3953eefcdc
# ╠═99420d2d-2eaf-4dd3-9dc7-ea08fa955be5
# ╠═368ce9a1-4ec8-4f9c-a148-68325a4b7194
# ╠═18ad0ba7-0d61-49c2-b9eb-9ef2a739e40d
# ╠═903130b1-980b-4259-868a-66f207d48e41
# ╠═47c10091-3b39-4f6a-bdde-0b1bbcf10e30
# ╠═0d13c2e4-59c7-43fb-a918-0715f232a406
# ╠═4d292003-38af-4060-bd5c-6a6a1e2244f1
# ╠═cc9b6619-1f7d-462e-962d-5027eb8b5c3c
# ╠═823dd2f5-d3bd-4069-bbe1-ad838ddca4a8
# ╠═bd0905fe-00bf-4ad0-b0d1-5d5d9aba04cf
# ╠═2a06b4e8-b3e8-436f-848b-bde664175505
# ╠═4b9000b8-2210-48ae-afef-08055fa1851d
# ╠═374a89d2-708f-4d6b-9daa-49363a68f218
# ╠═5a9d5d42-ddda-4b2e-aeae-020a39ae7c35
# ╠═14a0a868-c6ca-40a6-95fd-4e2052936ca3
# ╠═a5aadcc0-8b47-4dde-8373-800d4583e369
# ╠═3d058953-c4a8-453b-b0ca-7ad8fa930c48
# ╠═0ad36213-5639-443e-81d3-d8904f507639
# ╠═e501ad52-3ed5-4151-8dc1-7a362158e8c3
# ╠═3b2038e7-87a7-4bd7-b823-f230058c95b9
# ╠═eb74c6e9-1db0-4dbc-8222-47b5fc8744db
# ╠═96f6d024-5c55-4b3b-adbf-4cc8eb996adb
# ╠═e0c2c9a3-84be-4857-8e2a-022e29a1ae77
# ╠═c4b7d46c-411c-462b-b326-d2816c61c8c9
# ╠═4b3a570c-ba16-4bda-b8d8-81ade4b5e3e7
# ╠═cc4fcea9-eefe-4f5c-95bc-a1eac6873767
# ╟─d25871e5-5e97-4b56-9c1c-0eac45d778d2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
