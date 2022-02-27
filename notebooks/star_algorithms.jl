### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 223497c9-ae2d-4b10-992e-8f387b34b094
using LinearAlgebra

# ╔═╡ 6f5bd250-9b64-423a-83d9-843a1e0f8ec1
using Random

# ╔═╡ ccbb6059-4696-4648-9f1e-a696990ca816
using LinearMaps

# ╔═╡ 344496aa-e778-4256-999f-97231aff8713
using IterativeSolvers

# ╔═╡ 712fa486-97a9-11ec-0362-573810d1c4b9
md"""
# Star product algorithms

This is a cross-comparison and validation of computing the star product in various manners using the associativity of the star product
"""

# ╔═╡ 9904f6f2-0ee2-4548-ac95-efca3b423f94
"Calculate the scattering matrix star product"
function ⋆ₛ(A::AbstractMatrix, B::AbstractMatrix)
    d = Int(size(A, 1) // 2)
    @assert size(A) == size(B) == reverse(size(A)) "maps must be endomorphisms"
    invI_A₂₂B₁₁ = inv(I - A[(d+1):2d, (d+1):2d] * B[1:d, 1:d])
    invI_B₁₁A₂₂ = inv(I - B[1:d, 1:d] * A[(d+1):2d, (d+1):2d])
    [
        A[1:d, 1:d] + A[1:d, (d+1):2d] * invI_B₁₁A₂₂ * B[1:d, 1:d] * A[(d+1):2d, 1:d]        A[1:d, (d+1):2d] * invI_B₁₁A₂₂ * B[1:d, (d+1):2d];
        B[(d+1):2d, 1:d] * invI_A₂₂B₁₁ * A[(d+1):2d, 1:d]        B[(d+1):2d, (d+1):2d] + B[(d+1):2d, 1:d] * invI_A₂₂B₁₁ * A[(d+1):2d, (d+1):2d] * B[1:d, (d+1):2d];
    ]
end

# ╔═╡ 190d702e-6b1b-4bc9-830b-b502f4907208
begin
	nesting_level(x::Tuple{<:Tuple,<:Tuple}) = max(nesting_level.(x)...) + 1
	nesting_level(x::Tuple{Int,<:Tuple}) = nesting_level(x[2]) + 1
	nesting_level(x::Tuple{<:Tuple,Int}) = nesting_level(x[1]) + 1
	nesting_level(::Tuple{Int,Int}) = 1
end

# ╔═╡ 1f4f1a45-5c85-4f1d-abb1-f8acba8bb647
md"""
## Verification

Here we just compare that various nesting algorithms give the same output despite different nesting patterns and depths
"""

# ╔═╡ 58db25d7-8b8e-43c6-8e3f-3051b69b7232
min_nesting(n::Int) = 2^Int(ceil(log2(n))-1)

# ╔═╡ e25fbb87-b456-4406-a66c-a6d3903c4026
med_nesting(n::Int) = Int(floor(n/2))

# ╔═╡ 9be451cb-a23a-4098-b194-4303ebc24e46
max_nesting(::Int) = 1

# ╔═╡ a134ac2c-6c82-4a73-8ec6-350e48e69d5d
function get_As(l; n=20, seed=222)
	Random.seed!(seed)
	[rand(n,n) for _ in 1:l]
end

# ╔═╡ 801d1e4f-6043-4a23-b9b1-712348e44f51
As = get_As(3);

# ╔═╡ 1780fd0c-bd0f-4868-9d7b-81bd7f57e6b1
md"""
## MVP with GMRES

In this section we repeat the previous calculations but only using the star product as a matrix-vector product (MVP)
"""

# ╔═╡ dc60d700-064d-4329-8e8c-ec8fd6bc9fcb
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
        y₂ = gmres(I - A₂₂ * B₁₁, y[(d+1):2d])
        z = B * vcat(fill(zero(eltype(x)), d), x₂)
        z₁ = gmres(I - B₁₁ * A₂₂, z[1:d])
        z₂ = @view z[(d+1):2d]
        vcat(
            y₁ + A₁₂ * (B₁₁ * y₂ + z₁),
            B₂₁ * (y₂ + A₂₂ * z₁) + z₂
        )
    end
end

# ╔═╡ a9bfe3bf-3c74-42bc-bdda-b5073a62c6c6
"""
A recursive function to compute multiple star products. Accepts `f`, a function that calculates at which index in the chain to split the chain by performing a star product. This determines the degree of nesting of the products, and accepts `As`, a vector of matrices on which to take the star product.
"""
function star_products(f::Function, As::Vector)
	if length(As) == 1
		return (nest=1, mat=As[1])
	elseif length(As) == 2
		return (nest=(1, 1), mat=As[1] ⋆ₛ As[2])
	else
		n::Int = f(length(As))
		l = star_products(f, As[1:n])
		r = star_products(f, As[(n+1):end])
		return (nest=(l.nest, r.nest), mat=l.mat ⋆ₛ r.mat)
	end
end

# ╔═╡ 472464af-0c32-4eab-a51c-5331e4b32068
min_res = star_products(min_nesting, As);

# ╔═╡ b4119cc6-958c-4054-832c-b3f589d78f29
(max_depth=nesting_level(min_res.nest), nest=min_res.nest)

# ╔═╡ c5e19959-28cf-4eb4-b8bf-3955e47f893b
med_res = star_products(med_nesting, As);

# ╔═╡ bf879f47-53d9-4674-88d4-99e9cb05bdaa
(max_depth=nesting_level(med_res.nest), nest=med_res.nest)

# ╔═╡ 1c952f57-b642-430a-b315-424535da2d30
max_res = star_products(max_nesting, As);

# ╔═╡ d5b479f9-866b-4084-99a9-5c6d5a5cea09
(max_depth=nesting_level(max_res.nest), nest=max_res.nest)

# ╔═╡ 9a0d6b8a-8438-409c-8424-8eea7252cd4d
min_res.mat ≈ med_res.mat ≈ max_res.mat

# ╔═╡ fbfb4718-28fc-4de6-bdc7-e65219f06182
Bs = LinearMap.(As);

# ╔═╡ 1d843920-75a7-4015-824c-bd9cdb0cd4bb
free_min_res = star_products(min_nesting, Bs)

# ╔═╡ b4db11ea-58dc-4b5c-9f2a-c1287d554e1c
free_med_res = star_products(med_nesting, Bs)

# ╔═╡ d9a8fff6-dccc-4932-b672-a32ad388ed02
free_max_res = star_products(max_nesting, Bs)

# ╔═╡ 2f22ffa6-6881-4bbd-80eb-06fea7aa5be6
begin
	Random.seed!(111)
	v = rand(20)
	(free_min_res.mat * v) ≈ (free_med_res.mat * v) ≈ (free_max_res.mat * v) ≈ (min_res.mat * v) ≈ (med_res.mat * v) ≈ (max_res.mat * v)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearMaps = "7a12625a-238d-50fd-b39a-03d52299707e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
IterativeSolvers = "~0.9.2"
LinearMaps = "~3.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearMaps]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "dbb14c604fc47aa4f2e19d0ebb7b6416f3cfa5f5"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.5.1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╟─712fa486-97a9-11ec-0362-573810d1c4b9
# ╠═223497c9-ae2d-4b10-992e-8f387b34b094
# ╠═6f5bd250-9b64-423a-83d9-843a1e0f8ec1
# ╠═9904f6f2-0ee2-4548-ac95-efca3b423f94
# ╠═a9bfe3bf-3c74-42bc-bdda-b5073a62c6c6
# ╠═190d702e-6b1b-4bc9-830b-b502f4907208
# ╟─1f4f1a45-5c85-4f1d-abb1-f8acba8bb647
# ╠═58db25d7-8b8e-43c6-8e3f-3051b69b7232
# ╠═e25fbb87-b456-4406-a66c-a6d3903c4026
# ╠═9be451cb-a23a-4098-b194-4303ebc24e46
# ╠═a134ac2c-6c82-4a73-8ec6-350e48e69d5d
# ╠═801d1e4f-6043-4a23-b9b1-712348e44f51
# ╠═472464af-0c32-4eab-a51c-5331e4b32068
# ╠═b4119cc6-958c-4054-832c-b3f589d78f29
# ╠═c5e19959-28cf-4eb4-b8bf-3955e47f893b
# ╠═bf879f47-53d9-4674-88d4-99e9cb05bdaa
# ╠═1c952f57-b642-430a-b315-424535da2d30
# ╠═d5b479f9-866b-4084-99a9-5c6d5a5cea09
# ╠═9a0d6b8a-8438-409c-8424-8eea7252cd4d
# ╟─1780fd0c-bd0f-4868-9d7b-81bd7f57e6b1
# ╠═ccbb6059-4696-4648-9f1e-a696990ca816
# ╠═344496aa-e778-4256-999f-97231aff8713
# ╠═dc60d700-064d-4329-8e8c-ec8fd6bc9fcb
# ╠═fbfb4718-28fc-4de6-bdc7-e65219f06182
# ╠═1d843920-75a7-4015-824c-bd9cdb0cd4bb
# ╠═b4db11ea-58dc-4b5c-9f2a-c1287d554e1c
# ╠═d9a8fff6-dccc-4932-b672-a32ad388ed02
# ╠═2f22ffa6-6881-4bbd-80eb-06fea7aa5be6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
