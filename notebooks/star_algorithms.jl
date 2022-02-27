### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 223497c9-ae2d-4b10-992e-8f387b34b094
using LinearAlgebra

# ╔═╡ 6f5bd250-9b64-423a-83d9-843a1e0f8ec1
using Random

# ╔═╡ 712fa486-97a9-11ec-0362-573810d1c4b9
md"""
# Star product algorithms

This is a cross-comparison and validation of computing the star product in various manners using the associativity of the star product
"""

# ╔═╡ 9904f6f2-0ee2-4548-ac95-efca3b423f94
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

# ╔═╡ 714101fb-4c4b-424f-b2f6-132ba7b808e1
function max_nested_alg(As)
	Atot = As[1]
	for A in As[2:end]
		Atot = Atot ⋆ₛ A
	end
	return Atot
end

# ╔═╡ 35c033c6-7d9f-4fbb-a245-89783078b823
function min_nested_alg(As::Vector{<:Matrix})
	if length(As) == 1
		return As[1]
	elseif length(As) == 2
		return As[1] ⋆ₛ As[2]
	else
		l = round(Int, length(As)/2)
		return min_nested_alg(As[1:l]) ⋆ₛ min_nested_alg(As[(l+1):end])
	end
end

# ╔═╡ a134ac2c-6c82-4a73-8ec6-350e48e69d5d
function get_As(l; n=20, seed=222)
	Random.seed!(seed)
	[rand(n,n) for _ in 1:l]
end

# ╔═╡ 801d1e4f-6043-4a23-b9b1-712348e44f51
As = get_As(4)

# ╔═╡ 8b91f177-51b1-40ea-86d0-aaf81d5520ad
max_res = max_nested_alg(As)

# ╔═╡ c8138534-c46b-4cbb-be32-2573be19cd3c
min_res = min_nested_alg(As)

# ╔═╡ 9a0d6b8a-8438-409c-8424-8eea7252cd4d
min_res ≈ max_res

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
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

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╠═712fa486-97a9-11ec-0362-573810d1c4b9
# ╠═223497c9-ae2d-4b10-992e-8f387b34b094
# ╠═6f5bd250-9b64-423a-83d9-843a1e0f8ec1
# ╠═9904f6f2-0ee2-4548-ac95-efca3b423f94
# ╠═714101fb-4c4b-424f-b2f6-132ba7b808e1
# ╠═35c033c6-7d9f-4fbb-a245-89783078b823
# ╠═a134ac2c-6c82-4a73-8ec6-350e48e69d5d
# ╠═801d1e4f-6043-4a23-b9b1-712348e44f51
# ╠═8b91f177-51b1-40ea-86d0-aaf81d5520ad
# ╠═c8138534-c46b-4cbb-be32-2573be19cd3c
# ╠═9a0d6b8a-8438-409c-8424-8eea7252cd4d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
