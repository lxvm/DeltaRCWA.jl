### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ e3e0b674-1dde-4bd5-a30b-813e542e01e5
using Symbolics

# ╔═╡ ed5a5cf8-9ae9-49b2-b4e6-c9e837e093cc
using LinearAlgebra

# ╔═╡ 5de05616-445d-4395-a1ea-fc0220a6d5f5
md"
# Symbolics.jl derivations of scattering matrices from Generalized Sheet Transition Conditions (GSTCs)

I did this by hand, but it's very complex and I would like to verify or at least
be able to reproduce my steps programatically
"

# ╔═╡ 75bd4b62-6244-4217-b38c-b171cccef531
md"""
## Extended functionality

``\newcommand{\bm}{\boldsymbol}``
"""

# ╔═╡ bcb4b693-b036-4226-bb44-04650b5e0f82
Base.:*(D::Differential, N::Num) = D(N)

# ╔═╡ ba030d89-49d3-43bb-8cde-622aed27303a
Base.:*(a::Number, x::Equation) = a*x.lhs ~ a*x.rhs

# ╔═╡ 13d67561-3599-4c53-89ff-6c241fd523be
Base.:-(x::Equation) = -1x

# ╔═╡ 4dfa8443-ef39-430e-b371-8235f714fa12
Base.:+(x::Equation, y::Equation) = (x.lhs + y.lhs) ~ (x.rhs + y.rhs)

# ╔═╡ 3de291b2-3935-4752-a132-80dd9c57cfa4
Base.zero(::Equation) = Num(0) ~ Num(0)

# ╔═╡ e2e236c3-b2ff-4386-827d-74d0786eacaf
Symbolics.SymbolicUtils.expand(x::Equation) = expand(x.lhs) ~ expand(x.rhs)

# ╔═╡ c9294f81-4b01-4d04-aa91-1218531b4d78
Symbolics.SymbolicUtils.simplify_fractions(x::Equation; kw...) = simplify_fractions(x.lhs; kw...) ~ simplify_fractions(x.rhs; kw...)

# ╔═╡ e66bd6f7-9586-40aa-a47e-b41b09a0d84f
latexify_md(args...; kwargs...) = Markdown.LaTeX(repr(MIME"text/latex"(), Symbolics.latexify(args...; kwargs...)))

# ╔═╡ 82ef72e5-edfc-4b5b-9ac9-6b4a7de16d01
my_isequal(x::Equation, y::Equation) = isequal(expand.((x.rhs-x.lhs, y.rhs-y.lhs))...)

# ╔═╡ 989216eb-919f-4635-bcd1-6837ab73efc4
md"""
## Stating the curl equations

Faraday's law and Ampère's law (respectively) say
```math
\begin{align}
-\mu \partial_t \bm{H} &= \nabla \times \bm{E}
\\
\epsilon \partial_t \bm{E} &= \nabla \times \bm{H}
\end{align}
```
"""

# ╔═╡ c3557c41-d121-4724-936f-fdb4dacea919
@variables H[1:3] E[1:3] r[1:3] t ϵ μ # r[1] = x; r[2] = y; r[3] = z

# ╔═╡ 5f5a4a19-2617-43aa-a3ff-a6e141741281
H⃗, E⃗, r⃗ = Symbolics.scalarize.((H, E, r))

# ╔═╡ 690079b0-d473-4056-8cd5-2c2650ecffc3
∂ₜ = Differential(t)

# ╔═╡ 5cfccc77-bb63-41ef-bdb8-361f2b8b61c4
∇ = Differential.(r⃗)

# ╔═╡ f3248de4-aaee-44b6-beba-9abab0f19849
Ampère = ϵ*∂ₜ.(E⃗) .~ ∇ × H⃗

# ╔═╡ 2ac56fe2-7b35-42c5-a521-37f3d689f576
Faraday = -μ*∂ₜ.(H⃗) .~ ∇ × E⃗

# ╔═╡ afa85ef1-6b08-4810-8a51-e91b09e1f0a4
md"""
## The curl equations in the planewave basis

Consider a plane wave at frequency $\omega$ with wavevector ``\bm{k}``:
```math
\exp(i(\bm{k} \cdot \bm{r} - \omega t))
```
By considering solutions in the planewave basis, we will find that there is a redundant field component we can eliminate from the curl equations to obtain an expression for the curl equations in terms of the tangential fields and the wave vectors.
The solutions from Ampère's law and Faraday's law should be equivalent assuming the dispersion relation (described later) is satisfied.
"""

# ╔═╡ 61c9e61e-dd23-47f2-92cd-61db4effa8ec
@variables k[1:3] ω i

# ╔═╡ ee6c745b-32fa-4a44-8de1-f1c51405c216
k⃗ = Symbolics.scalarize(k)

# ╔═╡ 60465ae6-5713-421b-b23b-da3f16670d91
pw = exp(i*(k⃗⋅r⃗ - ω*t))

# ╔═╡ f0f33aea-7296-4e49-a031-793c5047d50f
use_pw_basis(x) = substitute.(expand_derivatives.(substitute.(
	x, Ref(Dict([(el, el*pw) for el in vcat(H⃗, E⃗)]))
)), Ref(Dict(pw => 1, i => 1)))

# ╔═╡ 3bb1486a-e748-4e38-bcb2-6755b2871696
Ampère_pw, Faraday_pw = use_pw_basis.((Ampère, Faraday))

# ╔═╡ 5216e525-a5b7-4d4d-b75b-c30a76e77939
Ez_rule = Dict(E⃗[3] => Symbolics.solve_for(Ampère_pw[3], E⃗[3]));

# ╔═╡ 62dd1407-a08e-497d-9739-d513eae93224
Faraday_solution = simplify_fractions.(expand.(substitute.(-ω*ϵ*Faraday_pw[1:2], Ref(Ez_rule))))

# ╔═╡ 46e529b9-9089-46c2-bfc4-1f859b9daea7
Hz_rule = Dict(H⃗[3] => Symbolics.solve_for(Faraday_pw[3], H⃗[3]));

# ╔═╡ 38df7d76-8d4d-4bc5-b626-11e4dbc67388
Ampère_solution = simplify_fractions.(expand.(substitute.(ω*μ*Ampère_pw[1:2], Ref(Hz_rule))))

# ╔═╡ 20e938a0-4f0b-4229-bade-e6317517f867
md"""
We can verify that this gives the expected results for Faraday's Law and Ampère's Law (respectively) of
```math
\begin{align}
(\bm{R K R}^{-1} - \omega^2 \epsilon \mu \bm I) \tilde{\bm H}_\parallel
&= \omega \epsilon \bm{R}^{-1} k_z \tilde{\bm E}_\parallel
\\
(\bm{R K R}^{-1} - \omega^2 \epsilon \mu \bm I) \tilde{\bm E}_\parallel
&= \omega \mu \bm R k_z \tilde{\bm H}_\parallel
\end{align}
```
where
```math
\bm{R} =
\begin{pmatrix}
0 & -1
\\
1 & 0
\end{pmatrix}
,
\bm K =
\begin{pmatrix}
k_x k_x & k_x k_y
\\
k_y k_x & k_y k_y
\end{pmatrix}
.
```
"""

# ╔═╡ 67485c4e-78fc-4849-991d-ea410c2894a4
K = k⃗[1:2]*transpose(k⃗[1:2]);

# ╔═╡ b0ebb818-4fc0-4e99-9bd8-942eea32411d
R = [
	0 -1;
	1 0;
];

# ╔═╡ a3470ead-59ba-47dd-a38a-db416726cf68
K_term = K - ω^2*ϵ*μ*I;

# ╔═╡ 20ea85fb-9171-449f-a598-85f95c21830b
my_Faraday = R * K_term * R' * H⃗[1:2] .~ ω*ϵ*k[3]*R'*E⃗[1:2]

# ╔═╡ 37e8ee6e-c762-47e0-8867-2661e7c6ef35
my_Ampère = R * K_term * R' * E⃗[1:2] .~ ω*μ*k[3]*R*H⃗[1:2]

# ╔═╡ 03c7d185-b0de-417e-96a9-495f45e6058f
all(my_isequal.(my_Ampère, Ampère_solution)) && all(my_isequal.(my_Faraday, Faraday_solution))

# ╔═╡ 4b9b68bf-99aa-4bcc-b1f2-9cd626b01bd2
md"""
## Stating the GSTCs

```math
\begin{align}
\bm{\rho}^e \bm{R} (\bm{H}_\parallel^2 - \bm{H}_\parallel^1) &= (\bm{E}_\parallel^2 + \bm{E}_\parallel^1)/2
\\
\bm{R}^{-1} (\bm{E}_\parallel^2 - \bm{E}_\parallel^1) &= \bm{\sigma}^m (\bm{H}_\parallel^2 + \bm{H}_\parallel^1)/2
\end{align}
```
"""

# ╔═╡ 8f2e2e12-3585-4159-b2c6-afa39ebdb646
@variables ρ[1:2, 1:2], σ[1:2, 1:2]

# ╔═╡ eebb88f9-68e4-49e9-90e7-942bc88a42a7
ρᵉ, σᵐ = Symbolics.scalarize.([ρ, σ])

# ╔═╡ 3e46c6ab-860b-4b3e-829f-3812fe424958
@variables H¹[1:2] E¹[1:2] H²[1:2] E²[1:2]

# ╔═╡ 23c72eb8-da61-44d2-8a80-b31dab90ffc8
H⃗¹, E⃗¹, H⃗², E⃗² = Symbolics.scalarize.([H¹, E¹, H², E²])

# ╔═╡ 40a843f5-6a95-4617-83e4-527479b7617c
GSTCᵉ = ρᵉ * R * (H⃗² - H⃗¹) .~ (E⃗² + E⃗¹)/2

# ╔═╡ a20f965d-e5e8-47e7-ae1d-592fc047853f
GSTCᵐ = R' * (E⃗² - E⃗¹) .~ σᵐ * (H⃗² + H⃗¹)/2

# ╔═╡ 0da14cff-952e-4c79-8961-e7d4ea505788
md"""
## Direction of propagation of planewave modes

In the planewave basis, the wave equations have the dispersion relation
```math
\bm{k} \cdot \bm{k} = \omega^2 \epsilon \mu
.
```
Given the ``k_x, k_y`` components (which in DeltaRCWA are known from Bloch's theorem), this relation gives two solutions for the the ``k_z`` component.
It could be positive (``+``, propagating forward) or negative (``-``, propagating backward).
Let ``k_z`` is the solution to a quadratic equation and we have to
regulate the sign.
Hence each mode thus far actually 
We make the substitutions
```math
\begin{align}
\bm{H}^{i}_\parallel &\to \bm{H}^{(i,+)}_\parallel + \bm{H}^{(i,-)}_\parallel
\\
\bm{E}^{i}_\parallel &\to \bm{E}^{(i,+)}_\parallel + \bm{E}^{(i,-)}_\parallel
\end{align}
```
These signs give us the propagating directions, and thus the incident and scattered components.
"""

# ╔═╡ d36f5186-74ba-4e84-b6ce-f5da848d9290
@variables H¹⁺[1:2] H²⁺[1:2] H¹⁻[1:2] H²⁻[1:2]

# ╔═╡ efa74999-9bad-4ba6-a2fc-7ed3db43c769
@variables E¹⁺[1:2] E²⁺[1:2] E¹⁻[1:2] E²⁻[1:2]

# ╔═╡ 31814210-30ce-4a5d-be4c-51906c664cec
H⃗¹⁺, H⃗²⁺, H⃗¹⁻, H⃗²⁻ = Symbolics.scalarize.([H¹⁺, H²⁺, H¹⁻, H²⁻])

# ╔═╡ f82b6bf9-e993-4e0f-8630-94306df85ea8
E⃗¹⁺, E⃗²⁺, E⃗¹⁻, E⃗²⁻ = Symbolics.scalarize.([E¹⁺, E²⁺, E¹⁻, E²⁻])

# ╔═╡ c3471ccf-4fbb-4c6f-94b5-4af9296428af
propagation_rules = Dict(Iterators.flatten([[
	H⃗¹[j] => H⃗¹⁺[j] + H⃗¹⁻[j],
	H⃗²[j] => H⃗²⁺[j] + H⃗²⁻[j],
	E⃗¹[j] => E⃗¹⁺[j] + E⃗¹⁻[j],
	E⃗²[j] => E⃗²⁺[j] + E⃗²⁻[j],
] for j in 1:2]));

# ╔═╡ e8bfcd14-c2f1-42ea-aa12-35a02618c54d
GSTC_propagating = substitute.(vcat(GSTCᵉ, GSTCᵐ), Ref(propagation_rules))

# ╔═╡ 42200ab4-862c-4e9f-be25-e6b70ba20af2
md"""
## Eliminating the electric field

We will use our solutions to the planewave form of Faraday's law to eliminate the electric field from the GSTC, i.e. using

$(latexify_md.(Faraday_solution))

However, the ``k_3 = k_z`` terms originated from ``\partial_z`` in Maxwell's equations in the planewave basis, and so to respect the direction of propagation of the planewave we must make the following substitution on the terms which originated from the partial derivative
```math
k_z \bm{E}^{i}_\parallel \to k_z \bm{E}^{(i,+)}_\parallel - k_z \bm{E}^{(i,-)}_\parallel
```
Doing so will give us the final form of the equations that we will interpret as a scattering matrix.
We could also use Ampère's law to do the same.
"""

# ╔═╡ 5376e9e8-565b-4107-ae8b-02af4b17b56a
E⃗_Faraday = Symbolics.solve_for.(Faraday_solution, E⃗[[2, 1]])[[2, 1]]

# ╔═╡ 821ccdc6-e804-42f3-a1fb-aa384ec3d4c3
E⃗_Faraday_rules = Dict(Iterators.flatten([
	zip(e, substitute.(E⃗_Faraday, Ref(Dict(zip(H⃗, h)))))
	for (e, h) in [(E⃗¹⁺, H⃗¹⁺), (E⃗²⁺, H⃗²⁺), (E⃗¹⁻, H⃗¹⁻), (E⃗²⁻, H⃗²⁻)]
]));

# ╔═╡ 3e86eb1e-30e1-446f-b7f7-d1ce7c23dc69
E⃗_Ampère = simplify.(Symbolics.solve_for(Ampère_solution, E⃗[1:2]), expand=true)

# ╔═╡ 83dab0da-f374-4f2b-b0cf-bdb54ae60fb7
E⃗_Ampère_rules = Dict(Iterators.flatten([
	zip(e, substitute.(E⃗_Ampère, Ref(Dict(zip(H⃗, h)))))
	for (e, h) in [(E⃗¹⁺, H⃗¹⁺), (E⃗²⁺, H⃗²⁺), (E⃗¹⁻, H⃗¹⁻), (E⃗²⁻, H⃗²⁻)]
]));

# ╔═╡ 98145011-3c0a-4c9f-aa2f-88f10e670f9f
GSTC_smatrix = expand.(substitute.(k[3].*GSTC_propagating, Ref(E⃗_Faraday_rules)))

# ╔═╡ 789829ca-0cd5-46df-b3c4-467dffea0f85
md"""
## Verifying the scattering matrix solution

These final 4 equations contain all the variables and information we need to calculate
the 4 scattered components from the 4 incident components.
Theoretically, we could do
```julia
Symbolics.solve_for(GSTC_smatrix, [H⃗¹⁻..., H⃗²⁺...]))
```
but this is unnecessary because we will use a matrix solver to do this in DeltaRCWA.
All we need to do in this notebook is to check that the equations above match what I got for the scattering matrix from my derivation, which I concisely express as
```math
\begin{align}
&\begin{pmatrix}
-(\bm K - \omega^2 \epsilon \mu \bm I) \bm R^{-1} \tilde{\bm \rho}^e \bm R
+ \frac{1}{2} \omega \mu k_z \bm I
&
(\bm K - \omega^2 \epsilon \mu \bm I) \bm R^{-1} \tilde{\bm \rho}^e \bm R
- \frac{1}{2} \omega \mu k_z \bm I
\\
-(\bm K - \omega^2 \epsilon \mu \bm I) \tilde{\bm \sigma}^m
+ 2\omega \mu k_z \bm I
&
-(\bm K - \omega^2 \epsilon \mu \bm I) \tilde{\bm \sigma}^m
+ 2\omega \mu k_z \bm I
\end{pmatrix}
\begin{pmatrix}
\tilde{\bm H}^{(1,-)}_\parallel
\\
\tilde{\bm H}^{(2,+)}_\parallel
\end{pmatrix}
\nonumber \\ =
&\begin{pmatrix}
(\bm K - \omega^2 \epsilon \mu \bm I) \bm R^{-1} \tilde{\bm \rho}^e \bm R
+ \frac{1}{2} \omega \mu k_z \bm I
&
-(\bm K - \omega^2 \epsilon \mu \bm I) \bm R^{-1} \tilde{\bm \rho}^e \bm R
- \frac{1}{2} \omega \mu k_z \bm I
\\
(\bm K - \omega^2 \epsilon \mu \bm I) \tilde{\bm \sigma}^m
+ 2\omega \mu k_z \bm I
&
(\bm K - \omega^2 \epsilon \mu \bm I) \tilde{\bm \sigma}^m
+ 2\omega \mu k_z \bm I
\end{pmatrix}
\begin{pmatrix}
\tilde{\bm H}^{(1,+)}_\parallel
\\
\tilde{\bm H}^{(2,-)}_\parallel
\end{pmatrix}
\end{align}
```
"""

# ╔═╡ f5838af8-bf11-46e3-883e-e229b774d8bf
ρᵉ_term = K_term * R' * ρᵉ * R;

# ╔═╡ 21d4aaf0-a547-4b32-96fc-2d5fa341ec16
σᵐ_term = K_term * σᵐ;

# ╔═╡ 79db383f-9710-4628-8c0f-eb220f4b76aa
LHS = [
	-ρᵉ_term + ω*μ*k[3]/2*I ρᵉ_term - ω*μ*k[3]/2*I;
	-σᵐ_term + 2ω*μ*k[3]*I -σᵐ_term + 2ω*μ*k[3]*I
];

# ╔═╡ ba4a3d8c-c50a-4a2a-accf-d7bad62225e4
RHS = [
	ρᵉ_term + ω*μ*k[3]/2*I -ρᵉ_term - ω*μ*k[3]/2*I;
	σᵐ_term + 2ω*μ*k[3]*I σᵐ_term + 2ω*μ*k[3]*I
];

# ╔═╡ b0d8bace-1005-491e-90d1-51cbd0c99150
my_smatrix = LHS * vcat(H⃗¹⁻[1:2], H⃗²⁺[1:2]) .~ RHS * vcat(H⃗¹⁺[1:2], H⃗²⁻[1:2])

# ╔═╡ 31659fdb-c809-4374-a29e-f8d28ab80504
all(my_isequal.(my_smatrix, GSTC_smatrix))

# ╔═╡ fdaa3f48-4cf4-4502-864c-23fef05a6650
md"""
Unfortunately I don't think we can trust checks for symbolic equality because functions like `simplify_fractions` did not cancel all instances of `ωϵ/(ωϵ)`, because the terms may be on different sides of the equal sign, because the equations may differ by a multiplicative constant, and mainly because in my calculation I eliminated the electric field using the planewave form of Ampère's law whereas in this symbolic calculation I've used Faraday's law.
(The two can be shown to be the same.)
I leave it to the reader to verify that term by term `my_GSTC` and `H_GSTC` are the same, or since it is a linear function of 4 inputs and 4 outputs, to check that numerically the results from each of the expressions agree on a basis for an 8d space.

For example, to verify the expressions for `E⃗_Faraday` and `E⃗_Ampère` are the same, one should use the rules below (and others like it) to go between them, however this kind of identity is more easily shown by hand than by Symbolics.jl.
"""

# ╔═╡ 0870f6b3-4db0-4998-b3f9-1f697b261b4b
k_rules = Dict(
	sum(k⃗.^2) => ω^2*ϵ*μ,
	-k[2]^2 + ω^2*μ*ϵ => k[1]^2 + k[3]^2,
);

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Symbolics = "~4.0.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgCheck]]
git-tree-sha1 = "dedbbb2ddb876f899585c4ec4433265e3017215a"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.1.0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "e527b258413e0c6d4f66ade574744c94edef81f8"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.40"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "0ad226aa72d8671f20d0316e03028f0ba1624307"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.32"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

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
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DefineSingletons]]
git-tree-sha1 = "77b4ca280084423b728662fe040e5ff8819347c5"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.1"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "3287dacf67c3652d3fed09f4c12c187ae4dbb89a"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.4.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "72dcda9e19f88d09bf21b5f9507a0bb430bce2aa"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.24"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "1b4665a7e303eaa7e03542cfaef0730cb056cb00"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.3.21"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "9aad812fb7c4c038da7cab5a069f502e6e3ae030"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.1"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[InitialValues]]
git-tree-sha1 = "7f6a4508b4a6f46db5ccd9799a3fc71ef5cad6e6"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.2.11"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "f0c6489b12d28fb4c2103073ec7452f3423bd308"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.1"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

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
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[LabelledArrays]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "8f5fd068dfee92655b79e0859ecad8b492dfe8b1"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.6.5"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "6193c3815f13ba1b78a51ce391db8be016ae9214"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.4"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0d3b2feb3168e4deb78361d3b5bb5c2e51ea5271"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.2"

[[MicroCollections]]
deps = ["BangBang", "Setfield"]
git-tree-sha1 = "4f65bdbbe93475f6ff9ea6969b21532f88d359be"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "45c9940cec79dedcdccc73cc6dd09ea8b8ab142c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.3.18"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "8d9496b2339095901106961f44718920732616bb"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.22"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c8b8775b2f242c80ea85c83714c64ecfa3c53355"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.3"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
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
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "c944fa4adbb47be43376359811c0a14757bdc8a8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.20.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "ad2c7f08e332cc3bb05d33026b71fa0ef66c009a"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.19.4"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "def0718ddbabeb5476e51e5a43609bee889f285d"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.0"

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
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "e7bc80dc93f50857a5d1e3c8121495852f407e6a"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "eb35dcc66558b2dda84079b9a1be17557d32091a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.12"

[[StatsFuns]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "95072ef1a22b057b1e80f73c2a89ad238ae4cfff"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.12"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "5a25c36acf4f288c88b4a83e3b87f0adc5e1f9f3"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.18.0"

[[Symbolics]]
deps = ["ConstructionBase", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "671779c01f26efbfc671edb924e085793ce301e6"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.0.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TermInterface]]
git-tree-sha1 = "897e35234f810b443868eb53873dfebb83998a0a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.2"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "abcff3ac31c7894550566be533b512f8b059104f"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.8"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "7cb456f358e8f9d102a8b25e8dfedf58fa5689bc"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.13"

[[Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "bccb153150744d476a6a8d4facf5299325d5a442"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.67"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─5de05616-445d-4395-a1ea-fc0220a6d5f5
# ╠═e3e0b674-1dde-4bd5-a30b-813e542e01e5
# ╠═ed5a5cf8-9ae9-49b2-b4e6-c9e837e093cc
# ╟─75bd4b62-6244-4217-b38c-b171cccef531
# ╠═bcb4b693-b036-4226-bb44-04650b5e0f82
# ╠═ba030d89-49d3-43bb-8cde-622aed27303a
# ╠═13d67561-3599-4c53-89ff-6c241fd523be
# ╠═4dfa8443-ef39-430e-b371-8235f714fa12
# ╠═3de291b2-3935-4752-a132-80dd9c57cfa4
# ╠═e2e236c3-b2ff-4386-827d-74d0786eacaf
# ╠═c9294f81-4b01-4d04-aa91-1218531b4d78
# ╠═e66bd6f7-9586-40aa-a47e-b41b09a0d84f
# ╠═82ef72e5-edfc-4b5b-9ac9-6b4a7de16d01
# ╟─989216eb-919f-4635-bcd1-6837ab73efc4
# ╠═c3557c41-d121-4724-936f-fdb4dacea919
# ╠═5f5a4a19-2617-43aa-a3ff-a6e141741281
# ╠═690079b0-d473-4056-8cd5-2c2650ecffc3
# ╠═5cfccc77-bb63-41ef-bdb8-361f2b8b61c4
# ╠═f3248de4-aaee-44b6-beba-9abab0f19849
# ╠═2ac56fe2-7b35-42c5-a521-37f3d689f576
# ╟─afa85ef1-6b08-4810-8a51-e91b09e1f0a4
# ╠═61c9e61e-dd23-47f2-92cd-61db4effa8ec
# ╠═ee6c745b-32fa-4a44-8de1-f1c51405c216
# ╠═60465ae6-5713-421b-b23b-da3f16670d91
# ╠═f0f33aea-7296-4e49-a031-793c5047d50f
# ╠═3bb1486a-e748-4e38-bcb2-6755b2871696
# ╠═5216e525-a5b7-4d4d-b75b-c30a76e77939
# ╠═62dd1407-a08e-497d-9739-d513eae93224
# ╠═46e529b9-9089-46c2-bfc4-1f859b9daea7
# ╠═38df7d76-8d4d-4bc5-b626-11e4dbc67388
# ╟─20e938a0-4f0b-4229-bade-e6317517f867
# ╠═67485c4e-78fc-4849-991d-ea410c2894a4
# ╠═b0ebb818-4fc0-4e99-9bd8-942eea32411d
# ╠═a3470ead-59ba-47dd-a38a-db416726cf68
# ╠═20ea85fb-9171-449f-a598-85f95c21830b
# ╠═37e8ee6e-c762-47e0-8867-2661e7c6ef35
# ╠═03c7d185-b0de-417e-96a9-495f45e6058f
# ╟─4b9b68bf-99aa-4bcc-b1f2-9cd626b01bd2
# ╠═8f2e2e12-3585-4159-b2c6-afa39ebdb646
# ╠═eebb88f9-68e4-49e9-90e7-942bc88a42a7
# ╠═3e46c6ab-860b-4b3e-829f-3812fe424958
# ╠═23c72eb8-da61-44d2-8a80-b31dab90ffc8
# ╠═40a843f5-6a95-4617-83e4-527479b7617c
# ╠═a20f965d-e5e8-47e7-ae1d-592fc047853f
# ╟─0da14cff-952e-4c79-8961-e7d4ea505788
# ╠═d36f5186-74ba-4e84-b6ce-f5da848d9290
# ╠═efa74999-9bad-4ba6-a2fc-7ed3db43c769
# ╠═31814210-30ce-4a5d-be4c-51906c664cec
# ╠═f82b6bf9-e993-4e0f-8630-94306df85ea8
# ╠═c3471ccf-4fbb-4c6f-94b5-4af9296428af
# ╠═e8bfcd14-c2f1-42ea-aa12-35a02618c54d
# ╟─42200ab4-862c-4e9f-be25-e6b70ba20af2
# ╠═5376e9e8-565b-4107-ae8b-02af4b17b56a
# ╠═821ccdc6-e804-42f3-a1fb-aa384ec3d4c3
# ╠═3e86eb1e-30e1-446f-b7f7-d1ce7c23dc69
# ╠═83dab0da-f374-4f2b-b0cf-bdb54ae60fb7
# ╠═98145011-3c0a-4c9f-aa2f-88f10e670f9f
# ╟─789829ca-0cd5-46df-b3c4-467dffea0f85
# ╠═f5838af8-bf11-46e3-883e-e229b774d8bf
# ╠═21d4aaf0-a547-4b32-96fc-2d5fa341ec16
# ╠═79db383f-9710-4628-8c0f-eb220f4b76aa
# ╠═ba4a3d8c-c50a-4a2a-accf-d7bad62225e4
# ╠═b0d8bace-1005-491e-90d1-51cbd0c99150
# ╠═31659fdb-c809-4374-a29e-f8d28ab80504
# ╟─fdaa3f48-4cf4-4502-864c-23fef05a6650
# ╠═0870f6b3-4db0-4998-b3f9-1f697b261b4b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
