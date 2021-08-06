### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 80095175-f0ad-4645-9e49-b1e83cf18c73
using Symbolics

# ╔═╡ 5de05616-445d-4395-a1ea-fc0220a6d5f5
md"
# Symbolic derivations of Generalized Sheet Transition Conditions (GSTCs)

I did this by hand, but it's very complex and I would like to verify or at least
be able to reproduce my steps with instructions in code
"

# ╔═╡ 989216eb-919f-4635-bcd1-6837ab73efc4
md"
## Stating the curl equations

Ampère's law says

$$-\mu \partial_t H_x = \partial_y E_z - \partial_z E_y$$
$$-\mu \partial_t H_y = \partial_z E_x - \partial_x E_z$$
$$-\mu \partial_t H_z = \partial_x E_y - \partial_y E_z$$

and Faraday's law says

$$\epsilon \partial_t E_x = \partial_y H_z - \partial_z H_y$$
$$\epsilon \partial_t E_y = \partial_z H_x - \partial_x H_z$$
$$\epsilon \partial_t E_z = \partial_x H_y - \partial_y H_z$$
"

# ╔═╡ f44f8519-b794-4f49-a926-378a6c33fc3a
# This is an attempt to vectorize the calculation, but the symbolic array variable
# declarations cause Julia to crash

# @variables H[1:3] E[1:3] r[0:3] # r[0] = t; r[1] = x; r[2] = y; r[3] = z
# vcurl_eqns = [
# 	Equation(
# 		F[1] * Differential(r[0])(F[2][ε[1]]),
# 		Differential(r[ε[2]])(F[3][ε[3]]) - Differential(r[ε[3]])(F[3][ε[2]])
# 	)
# 	for F in [(-μ, H, E),  (ϵ, E, H)]
# 	for ε in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
# ]
# @variables K[1:3]

# ╔═╡ f5311f88-fc9b-4c65-bd7f-bc1da71f00de
@variables x y z t

# ╔═╡ 7540fc36-dde4-499d-a7f2-21886f21c270
@variables Hˣ Hʸ Hᶻ Eˣ Eʸ Eᶻ

# ╔═╡ 39d79f4e-0aea-4102-b7e5-7125489de3b0
begin
	Dˣ = Differential(x)
	Dʸ = Differential(y)
	Dᶻ = Differential(z)
	Dᵗ = Differential(t)
end

# ╔═╡ afa85ef1-6b08-4810-8a51-e91b09e1f0a4
md"
## Solving the curl equations

We can simplify them to obtain equations only in terms of the $y$ field components

Consider a plane wave at frequency $\omega$ with wavevector $k_x, k_y, k_z$:

$\exp(i(k_x x + k_y y + k_z z - \omega t))$
"

# ╔═╡ dc9ce90c-e40b-47a3-9e25-e377e9d7875e
@variables kˣ kʸ kᶻ ω ϵ μ i

# ╔═╡ d9cb3504-21e4-41d6-85bf-264af6a8b9ef
curl_eqns = [
	Equation(-μ * Dᵗ(Hˣ), Dʸ(Eᶻ) - Dᶻ(Eʸ)),
	Equation(-μ * Dᵗ(Hʸ), Dᶻ(Eˣ) - Dˣ(Eᶻ)),
	Equation(-μ * Dᵗ(Hᶻ), Dᶻ(Eʸ) - Dʸ(Eᶻ)),
	Equation( ϵ * Dᵗ(Eˣ), Dʸ(Hᶻ) - Dᶻ(Hʸ)),
	Equation( ϵ * Dᵗ(Eʸ), Dᶻ(Hˣ) - Dˣ(Hᶻ)),
	Equation( ϵ * Dᵗ(Eᶻ), Dᶻ(Hʸ) - Dʸ(Hᶻ)),
]

# ╔═╡ 52200690-0f20-4ec7-a17d-5dcfe67297dd
# Nonlinear equations don't yet work

# dispersion = Equation(kˣ^2 + kʸ^2 + kᶻ^2, ω^2 * ϵ * μ)
# Symbolics.solve_for(dispersion, kᶻ)
# sqrt(ω^2 * ϵ * μ - kˣ^2 - kʸ^2)

# ╔═╡ 670a4092-3c54-4bbb-89a3-5e94f4beb0e7
planewave = exp(i * (kˣ * x + kʸ * y + kᶻ * z - ω * t))

# ╔═╡ 91a5c8ee-87b1-4032-9d67-fa7536d703ef
md"
We will use the planewave form of the solution to solve for the field components.
"

# ╔═╡ 5c434558-1aca-4391-9788-5af771cdc9e4
@variables hˣ hʸ hᶻ eˣ eʸ eᶻ

# ╔═╡ a8a444d9-8f99-4a0d-83d1-a05b36d0eff9
algebraic_curl_eqns = @. substitute(expand_derivatives(substitute(curl_eqns, Dict(
	Hˣ => hˣ * planewave,
	Hʸ => hʸ * planewave,
	Hᶻ => hᶻ * planewave,
	Eˣ => eˣ * planewave,
	Eʸ => eʸ * planewave,
	Eᶻ => eᶻ * planewave
))), Dict(planewave => 1, i => 1))
# Note that simplify will not cancel all common factors (issue below)

# ╔═╡ 25a2aabe-0127-409e-8ada-f3e444a911d8
# Maybe raise an issue in Symbolics.jl? MWE:
# @variables x z
# simplify((exp(z) + 2 * exp(z))/exp(z)) # gives 3
# simplify((exp(z) + 2x * exp(z))/exp(z)) # should give 2x + 1

# ╔═╡ 90966461-af5f-4f0d-bfac-7adf8c4967a1
coefs = [hˣ hʸ hᶻ eˣ eʸ eᶻ]

# ╔═╡ 087c47e1-fc39-45f4-addf-6967f53fb5a1
curled_coefs = [Symbolics.solve_for(algebraic_curl_eqns[j], coefs[j])
	for j in eachindex(coefs)]

# ╔═╡ ca713018-31d4-4e8b-bf4a-35c91bb5b95a
begin
	x_ind = [1, 4]
	y_ind = [2, 5]
	z_ind = [3, 6]
end

# ╔═╡ f61a33c9-ff96-4693-ba0c-c0f908148d88
# solve for hᶻ, eᶻ in terms of hʸ, eʸ
z_coefs_sol = Dict([
	(coefs[j[1]],
		Symbolics.solve_for(substitute(
			algebraic_curl_eqns[j[1]],
			Dict(coefs[j[2]] => curled_coefs[j[2]])
		), coefs[j[1]])
	)
	for j in zip(z_ind, reverse(z_ind))
])

# ╔═╡ b68afd45-3073-4bc0-ba51-3ba8387b153f
# solve for hˣ, eˣ in terms of hʸ, eʸ
x_coefs_sol = Dict([
	(coefs[j[1]],
		Symbolics.solve_for(substitute(
			algebraic_curl_eqns[j[1]],
			z_coefs_sol
		), coefs[j[1]])
	)
	for j in zip(x_ind, reverse(z_ind))
])

# ╔═╡ 4b9b68bf-99aa-4bcc-b1f2-9cd626b01bd2
md"
## Stating the GSTCs

We can give the GSTCs as a set of 4 equations

$$-[[H_y]] = \frac{\sigma_e \cdot \hat{x}}{2} \{\!\{ E_x \}\!\}$$
$$[[H_x]] = \frac{\sigma_e \cdot \hat{y}}{2} \{\!\{ E_y \}\!\}$$
$$-[[E_y]] = -\frac{\sigma_m \cdot \hat{x}}{2} \{\!\{ H_x \}\!\}$$
$$[[E_x]] = -\frac{\sigma_m \cdot \hat{y}}{2} \{\!\{ H_y \}\!\}$$

where double brackets give the difference between fields on each side of the sheet
and double braces give the sum of fields on each side of the sheet
"

# ╔═╡ 6bc79da4-f7b3-4f3c-ae3d-017ddc05d60c
@variables Hˣ₁ Eˣ₁ Hˣ₂ Eˣ₂

# ╔═╡ 5773d655-d338-4549-a065-037b8d408a45
@variables Hʸ₁ Eʸ₁ Hʸ₂ Eʸ₂

# ╔═╡ 585ecd23-4359-44a9-9e38-9a51fd856707
@variables σₑˣˣ σₑˣʸ σₑʸˣ σₑʸʸ

# ╔═╡ c0cecaa6-d832-430f-8b7e-2c28b6a67870
@variables σₘˣˣ σₘˣʸ σₘʸˣ σₘʸʸ

# ╔═╡ 9ebcd410-30d0-4b28-b505-fb1593a4e2d5
begin
	GSTCs = [
		Equation(-(Hʸ₁ - Hʸ₂), (σₑˣˣ/2) * (Eˣ₁ + Eˣ₂) + (σₑˣʸ/2) * (Eʸ₁ + Eʸ₂)),
		Equation((Hˣ₁ - Hˣ₂), (σₑʸʸ/2) * (Eʸ₁ + Eʸ₂) + (σₑʸˣ/2) * (Eˣ₁ + Eˣ₂)),
		Equation(-(Eʸ₁ - Eʸ₂), -(σₘˣˣ/2) * (Hˣ₁ + Hˣ₂) - (σₘˣʸ/2) * (Hʸ₁ + Hʸ₂)),
		Equation((Eˣ₁ - Eˣ₂), -(σₘʸʸ/2) * (Hʸ₁ + Hʸ₂) - (σₘʸˣ/2) * (Hˣ₁ + Hˣ₂)),
	]
end

# ╔═╡ 388cfd86-8b2e-41bc-b218-bd34d672beee
md"
Taking our previous solutions for the $x$ components, we can express it all with $y$
"

# ╔═╡ 04ad757e-5230-4015-876a-c1e373b5cd8f
@variables hʸ₁ eʸ₁ hʸ₂ eʸ₂

# ╔═╡ 94976229-bf33-4089-939b-e2ba9497be70
GSTCs_in_y = @. substitute(GSTCs, Dict(
	Hʸ₁ => hʸ₁,
	Hʸ₂ => hʸ₂,
	Eʸ₁ => eʸ₁,
	Eʸ₂ => eʸ₂,
	Eˣ₁ => substitute(x_coefs_sol[eˣ], Dict(eʸ => eʸ₁, hʸ => hʸ₁)),
	Hˣ₁ => substitute(x_coefs_sol[hˣ], Dict(eʸ => eʸ₁, hʸ => hʸ₁)),
	Eˣ₂ => substitute(x_coefs_sol[eˣ], Dict(eʸ => eʸ₂, hʸ => hʸ₂)),
	Hˣ₂ => substitute(x_coefs_sol[hˣ], Dict(eʸ => eʸ₂, hʸ => hʸ₂)),
))

# ╔═╡ 0da14cff-952e-4c79-8961-e7d4ea505788
md"
## Setup propagation directions of electromagnetic modes

The $k_z$ component could be traveling in the forward (+) or backward (-) directions,
whereas the $k_x, k_y$ components are expected to already have a sign built into them.
Recall this is because $k_z$ is the solution to a quadratic equation and we have to
regulate the sign.

These signs give us the propagating directions, and thus the incident and scattered
components. The $k_z$ terms originated from $\partial_z$ in Maxwell's equations.
"

# ╔═╡ bfb5cc1a-3580-48e8-b115-c2ec66ffa720
@variables hʸ₁_fw eʸ₁_fw hʸ₂_fw eʸ₂_fw

# ╔═╡ 1b40ef00-b7b8-440e-abe8-bccf5cea067d
@variables hʸ₁_bk eʸ₁_bk hʸ₂_bk eʸ₂_bk

# ╔═╡ 56de4210-610a-4ed2-9229-751ec55c6610
tmp_GSTCs_in_y = @. substitute(GSTCs_in_y, Dict(
	hʸ₁ => hʸ₁_fw + hʸ₁_bk,
	eʸ₁ => eʸ₁_fw + eʸ₁_bk,
	hʸ₂ => hʸ₂_fw + hʸ₂_bk,
	eʸ₂ => eʸ₂_fw + eʸ₂_bk
))

# ╔═╡ e178a1ed-3f61-4b42-8df6-cea0b79843a2
md"
Now we have to substitute
$k_z = \pm \sqrt{\omega^2 \epsilon \mu - k_x^2 - k_y^2}$
depending on the direction of propagation
"

# ╔═╡ a4785b44-893e-43f2-9f2c-6cf9911ff73a
scattering_GSTCs_in_y = @. substitute(tmp_GSTCs_in_y, Dict(
	kᶻ * (hʸ₁_fw + hʸ₁_bk)  => (hʸ₁_fw - hʸ₁_bk) * sqrt(ω^2 * ϵ * μ - kˣ^2 - kʸ^2),
	kᶻ * (eʸ₁_fw + eʸ₁_bk)  => (eʸ₁_fw - eʸ₁_bk) * sqrt(ω^2 * ϵ * μ - kˣ^2 - kʸ^2),
	kᶻ * (hʸ₂_fw + hʸ₂_bk)  => (hʸ₂_fw - hʸ₂_bk) * sqrt(ω^2 * ϵ * μ - kˣ^2 - kʸ^2),
	kᶻ * (eʸ₂_fw + eʸ₂_bk)  => (eʸ₂_fw - eʸ₂_bk) * sqrt(ω^2 * ϵ * μ - kˣ^2 - kʸ^2),
))

# ╔═╡ a03e8a81-bdaa-4c91-87a9-eb07f8e58c71
md"
## Obtaining the solution

These final 4 equations contain all the variables and information we need to calculate
the 4 scattered components from the 4 incident components
"

# ╔═╡ 34351308-96e3-4688-b288-70dcca5625e2
scattered_comp = [hʸ₁_bk, eʸ₁_bk, hʸ₂_fw, eʸ₂_fw]

# ╔═╡ 1836d98b-df52-452d-970b-ac3c852b146e
incident_comp = [hʸ₁_fw, eʸ₁_fw, hʸ₂_bk, eʸ₂_bk]

# ╔═╡ 0f6e6fea-09bb-4b68-8ab1-65ae318af5ea
# Running this cell should be done with caution (I let it go 10 minutes, incompleted)

# Symbolics.solve_for(scattering_GSTCs_in_y, scattered_comp)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Symbolics = "~1.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "a71d224f61475b93c9e196e83c17c6ac4dedacfa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.18"

[[Artifacts]]
deps = ["Pkg"]
git-tree-sha1 = "c30985d8821e0cd73870b17b0ed0ce6dc44cb744"
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.3.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f53ca8d41e4753c41cdafa6ec5f7ce914b34be54"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.10.13"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8e695f735fca77e9708e795eda62afdb869cbb70"
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.3.4+0"

[[CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "85d2d9e2524da988bffaf2a381864e20d2dae08d"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.2.1"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3889f646423ce91dd1055a76317e9a1d3a23fff1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.11"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics", "Test"]
git-tree-sha1 = "6cdd99d0b7b555f96f7cb05aa82067ee79e7aef4"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.2"

[[DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "5e47c4d652ea67652b7c5945c79c46472397d47f"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.3.18"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "8041575f021cba5a099a456b4163c9a08b566a02"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "25b9cc23ba3303de0ad2eac03f840de9104c9253"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[LabelledArrays]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "5e38cfdd771c34821ade5515f782fe00865d60b3"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.6.2"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LibGit2]]
deps = ["Printf"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "4ea90bd5d3985ae1f9a908bd4500ae88921c5ce7"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "45c9940cec79dedcdccc73cc6dd09ea8b8ab142c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.3.18"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9db77584158d0ab52307f8c04f8e7c08ca76b5b3"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.3+4"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Pkg]]
deps = ["Dates", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "UUIDs"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "2a76e8f24c67f3ebecaccefa8d4abd27db828407"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.14.9"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "86c5647b565873641538d8f812c04e4c9dbeb370"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.6.1"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1b7bf41258f6c5c9c31df8c1ba34c1fc88674957"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.2.2+2"

[[RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "5975a4f824533fa4240f40d86f1060b9fc80d7cc"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "932aaae93e81686e6473f27b5a11c403576a2183"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.18.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "d5640fc570fb1b6c54512f0bd3853866bd298b3e"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.0"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a50550fa3164a8c46747e62063b4d774ac1bcf49"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.5.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "1b9a0f17ee0adde9e538227de093467348992397"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.7"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2f6792d523d7448bbe2fec99eca9218f06cc746d"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.8"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SymbolicUtils]]
deps = ["AbstractTrees", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs"]
git-tree-sha1 = "17fecd52e9a82ca52bfc9d28a2d31b33458a6ce0"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.13.1"

[[Symbolics]]
deps = ["ConstructionBase", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "e1e41ac302c7a28875c1ab9917499a44d037362b"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "1.4.0"

[[TOML]]
deps = ["Dates"]
git-tree-sha1 = "44aaac2d2aec4a850302f9aa69127c74f0c3787e"
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "8ed4a3ea724dac32670b062be3ef1c1de6773ae8"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.4.4"

[[Test]]
deps = ["Distributed", "InteractiveUtils", "Logging", "Random"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "209a8326c4f955e2442c07b56029e88bb48299c7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.12"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "9e7a1e8ca60b742e508a315c17eef5211e7fbfd7"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.1"
"""

# ╔═╡ Cell order:
# ╟─5de05616-445d-4395-a1ea-fc0220a6d5f5
# ╠═80095175-f0ad-4645-9e49-b1e83cf18c73
# ╟─989216eb-919f-4635-bcd1-6837ab73efc4
# ╠═f44f8519-b794-4f49-a926-378a6c33fc3a
# ╠═f5311f88-fc9b-4c65-bd7f-bc1da71f00de
# ╠═7540fc36-dde4-499d-a7f2-21886f21c270
# ╠═39d79f4e-0aea-4102-b7e5-7125489de3b0
# ╠═d9cb3504-21e4-41d6-85bf-264af6a8b9ef
# ╟─afa85ef1-6b08-4810-8a51-e91b09e1f0a4
# ╠═dc9ce90c-e40b-47a3-9e25-e377e9d7875e
# ╠═52200690-0f20-4ec7-a17d-5dcfe67297dd
# ╠═670a4092-3c54-4bbb-89a3-5e94f4beb0e7
# ╟─91a5c8ee-87b1-4032-9d67-fa7536d703ef
# ╠═5c434558-1aca-4391-9788-5af771cdc9e4
# ╠═a8a444d9-8f99-4a0d-83d1-a05b36d0eff9
# ╠═25a2aabe-0127-409e-8ada-f3e444a911d8
# ╠═90966461-af5f-4f0d-bfac-7adf8c4967a1
# ╠═087c47e1-fc39-45f4-addf-6967f53fb5a1
# ╠═ca713018-31d4-4e8b-bf4a-35c91bb5b95a
# ╠═f61a33c9-ff96-4693-ba0c-c0f908148d88
# ╠═b68afd45-3073-4bc0-ba51-3ba8387b153f
# ╟─4b9b68bf-99aa-4bcc-b1f2-9cd626b01bd2
# ╠═6bc79da4-f7b3-4f3c-ae3d-017ddc05d60c
# ╠═5773d655-d338-4549-a065-037b8d408a45
# ╠═585ecd23-4359-44a9-9e38-9a51fd856707
# ╠═c0cecaa6-d832-430f-8b7e-2c28b6a67870
# ╠═9ebcd410-30d0-4b28-b505-fb1593a4e2d5
# ╟─388cfd86-8b2e-41bc-b218-bd34d672beee
# ╠═04ad757e-5230-4015-876a-c1e373b5cd8f
# ╠═94976229-bf33-4089-939b-e2ba9497be70
# ╟─0da14cff-952e-4c79-8961-e7d4ea505788
# ╠═bfb5cc1a-3580-48e8-b115-c2ec66ffa720
# ╠═1b40ef00-b7b8-440e-abe8-bccf5cea067d
# ╠═56de4210-610a-4ed2-9229-751ec55c6610
# ╟─e178a1ed-3f61-4b42-8df6-cea0b79843a2
# ╠═a4785b44-893e-43f2-9f2c-6cf9911ff73a
# ╟─a03e8a81-bdaa-4c91-87a9-eb07f8e58c71
# ╠═34351308-96e3-4688-b288-70dcca5625e2
# ╠═1836d98b-df52-452d-970b-ac3c852b146e
# ╠═0f6e6fea-09bb-4b68-8ab1-65ae318af5ea
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
