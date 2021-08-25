### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 61f4e85b-3b4e-4906-85a6-fd4a08056532
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ f4807a09-87dd-49c0-ab3d-ca6cb6cd046d
using Revise

# ╔═╡ ae0ef922-4287-4fb8-b3bb-f1cdfd4f3d6b
using DeltaRCWA

# ╔═╡ 7a6a3ee1-4fc8-4b1f-ba8a-5312fb2b4a1a
using Plots

# ╔═╡ 46dda2a7-9ca8-4c2e-b00e-e6aa95761f65
md"
# 3D problem

This notebook compares the DeltaRCWA 3D solver against the 2D one
"

# ╔═╡ 3afdc5c3-83c9-4630-a494-7dba45b03a00
struct ComplexExpSheet{T, N} <: RCWASheet{T, N}
    θ::T
    θᵗ::T
    k::T
    d::T
    L::T
end

# ╔═╡ be89853e-634d-469a-af19-822441f48b29
begin
	### extend DeltaRCWA with methods for this parametrized metasurface
	DeltaRCWA.ρₑˣˣ(sheet::ComplexExpSheet, x⃗) = 1e-8 .- 0.5sin(sheet.θ).*(1 .- exp.([im*sheet.k*sheet.d*e⃗[1] for e⃗ in Iterators.product(x⃗...)]))
	DeltaRCWA.σₑˣˣ(sheet::ComplexExpSheet, x⃗) = 1 ./ ρₑˣˣ(sheet, x⃗)
	DeltaRCWA.σₘʸʸ(sheet::ComplexExpSheet, x⃗) = -2sin(sheet.θ).*(1 .+ exp.([im*sheet.k*sheet.d*e⃗[1] for e⃗ in Iterators.product(x⃗...)]))
end;

# ╔═╡ 000ef374-c347-47eb-bd60-9e119af6ce11
begin
	### constants
	k = 10.0   # wavenumber k²=ω²μ₀ϵ₀
	ω = k + 1e-5 # do not choose ω = k because in prob2D this creates a value of kz==0
	θ = -π/2.0 # incidence angle (measured with respect to the x axis)
	θᵗ = -π/3    # transmitted field angle
	d  = cos(θᵗ)-cos(θ)
	L = 2*(2*π)/(k*abs(d));  # Unit cell width
	Nxmodes = 100
	Nymodes = 10
end;

# ╔═╡ d1ababf4-748f-499e-bd2e-42aa56a1ed96
begin
	### Define 1D problem
	sheet1D = ComplexExpSheet{Float64, 1}(θ, θᵗ, k, d, L)
	dims1D = ((Nxmodes, sheet1D.L), )
	pol1D = TM()
	I₁ = [n==1 ? 1 : 0 for n in 1:Nxmodes]
	I₂ = zeros(Nxmodes)
	prob1D = DeltaRCWAProblem(sheet1D, dims1D, ω, pol1D, I₁, I₂)
end;

# ╔═╡ a5bdc94e-b32d-4509-ba7a-cb1b05f17cd8
sol1D = solve(prob1D);

# ╔═╡ 4da6a0a2-417a-4a53-bb81-9915b7b557ac
plot(sol1D)

# ╔═╡ 1d3b118a-a8b0-4dad-aecb-a75e9ea3161b
begin
	### Define equivalent 2D problem
	sheet2D = ComplexExpSheet{Float64, 2}(θ, θᵗ, k, d, L)
	dims2D = ((Nxmodes, sheet2D.L), (Nymodes, sheet2D.L))
	pol2D = Coupled()
	I₁ˣ = [m==n==1 ? 1 : 0 for n in 1:Nxmodes, m in 1:Nymodes]
	I₂ˣ = zeros(Nxmodes, Nymodes)
	I₁ʸ = zeros(Nxmodes, Nymodes)
	I₂ʸ = zeros(Nxmodes, Nymodes)
	prob2D = DeltaRCWAProblem(sheet2D, dims2D, ω, pol2D, hcat(I₁ˣ, I₁ʸ), hcat(I₂ˣ, I₂ʸ))
end;

# ╔═╡ 896f11c8-1365-4eee-9fd3-a74e5c6e22b0
### if true, then there will be a singular s-matrix
any(iszero.(prob2D.modes.kz))

# ╔═╡ bf6c7112-5959-40d9-910f-005b2b732cf9
sol2D = solve(prob2D);

# ╔═╡ 689923b4-dc84-4fcb-bd1e-2404d3b0d5a4
# S2D = smatrix(sheet2D, prob2D.modes, pol2D);

# ╔═╡ 1ef94741-bac9-4195-8664-ee7791eac0b7
# heatmap(abs.(S2D))

# ╔═╡ 877eca8b-9e16-4a7d-ba56-8d44569ffaa3
begin
	# extract 1D slice of 2D solution for comparison
	i = 1::Int
	@assert 1 <= i <= Nymodes
	I₁slice_Hx = sol2D.I₁[:, i]
	I₁slice_Hy = sol2D.I₁[:, Nymodes+i]
	I₂slice_Hx = sol2D.I₂[:, i]
	I₂slice_Hy = sol2D.I₂[:, Nymodes+i]	
	O₁slice_Hx = sol2D.O₁[:, i]
	O₁slice_Hy = sol2D.O₁[:, Nymodes+i]
	O₂slice_Hx = sol2D.O₂[:, i]
	O₂slice_Hy = sol2D.O₂[:, Nymodes+i]
	Hx_slice_sol2D = DeltaRCWA.DeltaRCWASolution(prob1D.modes, pol1D, I₁slice_Hx, I₂slice_Hx, O₁slice_Hx, O₂slice_Hx)
	Hy_slice_sol2D = DeltaRCWA.DeltaRCWASolution(prob1D.modes, pol1D, I₁slice_Hy, I₂slice_Hy, O₁slice_Hy, O₂slice_Hy)
end;

# ╔═╡ b3d4c56f-48a6-4875-ac26-4d052f5d8bdc
plot(Hx_slice_sol2D)

# ╔═╡ 592db5ab-2255-4458-92b4-625505bb6b9f
plot(Hy_slice_sol2D)

# ╔═╡ Cell order:
# ╟─46dda2a7-9ca8-4c2e-b00e-e6aa95761f65
# ╠═61f4e85b-3b4e-4906-85a6-fd4a08056532
# ╠═f4807a09-87dd-49c0-ab3d-ca6cb6cd046d
# ╠═ae0ef922-4287-4fb8-b3bb-f1cdfd4f3d6b
# ╠═7a6a3ee1-4fc8-4b1f-ba8a-5312fb2b4a1a
# ╠═3afdc5c3-83c9-4630-a494-7dba45b03a00
# ╠═be89853e-634d-469a-af19-822441f48b29
# ╠═000ef374-c347-47eb-bd60-9e119af6ce11
# ╠═d1ababf4-748f-499e-bd2e-42aa56a1ed96
# ╠═a5bdc94e-b32d-4509-ba7a-cb1b05f17cd8
# ╠═4da6a0a2-417a-4a53-bb81-9915b7b557ac
# ╠═1d3b118a-a8b0-4dad-aecb-a75e9ea3161b
# ╠═896f11c8-1365-4eee-9fd3-a74e5c6e22b0
# ╠═bf6c7112-5959-40d9-910f-005b2b732cf9
# ╠═689923b4-dc84-4fcb-bd1e-2404d3b0d5a4
# ╠═1ef94741-bac9-4195-8664-ee7791eac0b7
# ╠═877eca8b-9e16-4a7d-ba56-8d44569ffaa3
# ╠═b3d4c56f-48a6-4875-ac26-4d052f5d8bdc
# ╠═592db5ab-2255-4458-92b4-625505bb6b9f
