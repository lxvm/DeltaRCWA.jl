script_dir = ENV["HOME"] * "/.julia/dev/DeltaRCWA/scripts/"
import Pkg
Pkg.activate(script_dir)

using FFTW
using LinearAlgebra
using Plots
using JLD2

include("Carlos_method.jl")
include("Luke_method.jl")
# include("Lorenzo_method.jl")

θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction
θᵗ = π/3 # desired transmitted field angle
d = cos(θᵗ)-cos(θⁱ) # design parameter
k₀ = 10.0 # wavenumber that attains the metasurface design goals
λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals
L₀ = λ₀/abs(d) # period of metasurface

η(x) = sin(abs(θⁱ)).*(1+exp.(im*k₀*d*x))
μ(x) = sin(abs(θⁱ)).*(1-exp.(im*k₀*d*x))

# η(x) = 0
# μ(x) = 0
L=1.0
k=1.1k₀

α = k * cos(θⁱ)
β = k * sin(θⁱ)
i = 1 # == 1 for θⁱ = π/2
exponents = 5:7
I₁error = zeros(length(exponents))
I₂error = zeros(length(exponents))
O₁error = zeros(length(exponents))
O₂error = zeros(length(exponents))
for e in eachindex(exponents)
    local N = 2^exponents[e]
    println("starting BIE #$e")
    local BIE = compute_BIE_method(η, μ, L, N, k, i)
    println("starting RCWA #$e")
    local RCWA = compute_RCWA_method(η, μ, L, N, k, i)
    # println("starting DeltaRCWA")
    # ΔRCWA = compute_DeltaRCWA_method(η, μ, L, N, k, i)
    I₁error[e] = sqrt(sum(abs2.(real.(BIE.I₁ - RCWA.I₁)))) / N^2
    O₁error[e] = sqrt(sum(abs2.(real.(BIE.O₁ - RCWA.O₁)))) / N^2
    I₂error[e] = sqrt(sum(abs2.(real.(BIE.I₂ - RCWA.I₂)))) / N^2
    O₂error[e] = sqrt(sum(abs2.(real.(BIE.O₂ - RCWA.O₂)))) / N^2
end
display(plot(exponents, I₁error, title="L2 norm of BIE vs RCWA incident fields below", xguide="log2(N)", yguide="error"))
display(plot(exponents, O₁error, title="L2 norm of BIE vs RCWA scattered fields below", xguide="log2(N)", yguide="error"))
display(plot(exponents, I₂error, title="L2 norm of BIE vs RCWA incident fields above", xguide="log2(N)", yguide="error"))
display(plot(exponents, O₂error, title="L2 norm of BIE vs RCWA scattered fields above", xguide="log2(N)", yguide="error"))