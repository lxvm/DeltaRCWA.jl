### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ f4a38678-e0e4-43e2-ad72-f2fd0900fdad
using IterativeSolvers

# ╔═╡ 093c22ac-dae6-4dd6-9064-1718545c8abf
using LinearAlgebra

# ╔═╡ a8b5972b-c41a-4d6a-aa3f-0bc4f3c63a4b
using KrylovKit

# ╔═╡ 183045de-fa33-11eb-141d-57f63848a09b
md"
# Iterative star product

The purpose of this is to implement sequential Redheffer star products
([wiki](https://en.wikipedia.org/wiki/Redheffer_star_product))
with iterative solvers.
"

# ╔═╡ 73180984-9230-4657-8975-466e1f4faff9
n=1000;

# ╔═╡ 5a2d417c-a188-40f7-b366-4ddc3e0b19a6
y = rand(n)

# ╔═╡ 6ac8fdc3-4b0b-4d34-851e-4ac7f9d787d9
# A = rand(n, n);
A = I(n);

# ╔═╡ 36dab249-491d-4432-a945-3cd5799fc27d
# time to solve Ax = y

# ╔═╡ 9ffa65fe-75fc-48d4-bb8a-159b6f8c0ae2
invLA = inv(A) * y

# ╔═╡ 60e9457e-be02-4823-bcc1-2e08c9009268
invgmres = gmres(A, y)

# ╔═╡ dd6bc5a1-11c0-4d7a-ba6d-c8f457eb7fbb
invKit, info = linsolve(A, y) 

# ╔═╡ 4f2b489d-a9df-4b3e-bac0-3e7f37e3eb88
invLA ≈ invgmres ≈ invKit

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
IterativeSolvers = "~0.9.1"
KrylovKit = "~0.5.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1a8c6237e78b714e901e406c096fc8a65528af7d"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.1"

[[KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "0328ad9966ae29ccefb4e1b9bfd8c8867e4360df"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.3"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═183045de-fa33-11eb-141d-57f63848a09b
# ╠═f4a38678-e0e4-43e2-ad72-f2fd0900fdad
# ╠═093c22ac-dae6-4dd6-9064-1718545c8abf
# ╠═a8b5972b-c41a-4d6a-aa3f-0bc4f3c63a4b
# ╠═73180984-9230-4657-8975-466e1f4faff9
# ╠═5a2d417c-a188-40f7-b366-4ddc3e0b19a6
# ╠═6ac8fdc3-4b0b-4d34-851e-4ac7f9d787d9
# ╠═36dab249-491d-4432-a945-3cd5799fc27d
# ╠═9ffa65fe-75fc-48d4-bb8a-159b6f8c0ae2
# ╠═60e9457e-be02-4823-bcc1-2e08c9009268
# ╠═dd6bc5a1-11c0-4d7a-ba6d-c8f457eb7fbb
# ╠═4f2b489d-a9df-4b3e-bac0-3e7f37e3eb88
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
