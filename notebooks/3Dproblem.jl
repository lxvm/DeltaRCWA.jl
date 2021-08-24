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
	DeltaRCWA.ρₑˣˣ(sheet::ComplexExpSheet, x⃗) = 2 .- 2sin(sheet.θ).*(1 .- exp.([im*sheet.k*sheet.d*e⃗[1] for e⃗ in Iterators.product(x⃗...)]))
	DeltaRCWA.σₑˣˣ(sheet::ComplexExpSheet, x⃗) = 1 ./ ρₑˣˣ(sheet, x⃗)
	DeltaRCWA.σₘʸʸ(sheet::ComplexExpSheet, x⃗) = -2sin(sheet.θ).*(1 .+ exp.([im*sheet.k*sheet.d*e⃗[1] for e⃗ in Iterators.product(x⃗...)]))
end;

# ╔═╡ 000ef374-c347-47eb-bd60-9e119af6ce11
begin
	### constants
	k = 10.0   # wavenumber k²=ω²μ₀ϵ₀
	ω = k
	θ = -π/2.0 # incidence angle (measured with respect to the x axis)
	θᵗ = -π/3    # transmitted field angle
	d  = cos(θᵗ)-cos(θ)
	L = 2*(2*π)/(k*abs(d));  # Unit cell width
	Nmodes = 100
end;

# ╔═╡ d1ababf4-748f-499e-bd2e-42aa56a1ed96
begin
	### Define 1D problem
	sheet1D = ComplexExpSheet{Float64, 1}(θ, θᵗ, k, d, L)
	dims1D = ((Nmodes, sheet1D.L), )
	pol1D = TM()
	I₁ = [n==1 ? 1 : 0 for n in 1:Nmodes]
	I₂ = zeros(Nmodes)
	prob1D = DeltaRCWAProblem(sheet1D, dims1D, ω, pol1D, I₁, I₂)
end;

# ╔═╡ a5bdc94e-b32d-4509-ba7a-cb1b05f17cd8
sol1D = solve(prob1D)

# ╔═╡ 4da6a0a2-417a-4a53-bb81-9915b7b557ac
plot(sol1D)

# ╔═╡ 1d3b118a-a8b0-4dad-aecb-a75e9ea3161b
begin
	### Define equivalent 2D problem
	sheet2D = ComplexExpSheet{Float64, 2}(θ, θᵗ, k, d, L)
	dims2D = ((Nmodes, sheet2D.L), (Nmodes, sheet2D.L))
	pol2D = Coupled()
	I₁ˣ = [m==n==1 ? 1 : 0 for n in 1:Nmodes, m in 1:Nmodes]
	I₂ˣ = zeros(Nmodes, Nmodes)
	I₁ʸ = zeros(Nmodes, Nmodes)
	I₂ʸ = zeros(Nmodes, Nmodes)
	prob2D = DeltaRCWAProblem(sheet2D, dims2D, ω, pol2D, hcat(I₁ˣ, I₁ʸ), hcat(I₂ˣ, I₂ʸ))
end;

# ╔═╡ 145c7554-77bd-4f8c-86c1-230bbd616c26
solve(prob2D)

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
# ╠═145c7554-77bd-4f8c-86c1-230bbd616c26
