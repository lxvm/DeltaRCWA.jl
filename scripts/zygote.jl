using DeltaRCWA
using LinearAlgebra
using Zygote

struct ComplexExpSheet{T} <: Sheet
    θ::T # incidence angle
    k::T # wavenumber
    d::T # designer parameter
end

function DeltaRCWA.Zₑˣˣ(sheet::ComplexExpSheet, x⃗::Tuple)
    -0.5sin(sheet.θ)*(1-exp(1im*sheet.k*sheet.d*x⃗[1]))
end

function DeltaRCWA.Yₘʸʸ(sheet::ComplexExpSheet, x⃗::Tuple)
    -2sin(sheet.θ)*(1+exp(1im*sheet.k*sheet.d*x⃗[1]))
end

function cost(sol::DeltaRCWA.DeltaRCWASolution)
    norm(sol.O₂) + norm(sol.O₁)
end

function grad2d()
    θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction
    θᵗ = π/3 # desired transmitted field angle
    d = cos(θᵗ)-cos(θⁱ) # design parameter
    k₀ = 10.0 # wavenumber that attains the metasurface design goals
    λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals
    L₀ = λ₀/abs(d) # period of metasurface
    sheet = ComplexExpSheet(θⁱ, k₀, d) # sheet to scatter

    L=L₀ # length of domain (some integer multiple of period of metasurface)
    k=1.1k₀ # wavenumber of incident field
    n = 20 # number of modes to scatter
    dims = ((n, L),)
    I₁ = [i == 1 ? 1.0 : 0.0 for i in 1:n]
    I₂ = zeros(n)
    gradient(s -> cost(solve(DeltaRCWAProblem(s, dims, k, I₁, I₂))), sheet)
end # ((θ = -3.831612187550387e-8, k = 6.62037566523508e10, d = 1.324075133047017e12),)

function grad3d()
    θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction
    θᵗ = π/3 # desired transmitted field angle
    d = cos(θᵗ)-cos(θⁱ) # design parameter
    k₀ = 10.0 # wavenumber that attains the metasurface design goals
    λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals
    L₀ = λ₀/abs(d) # period of metasurface
    sheet = ComplexExpSheet(θⁱ, k₀, d) # sheet to scatter

    L=L₀ # length of domain (some integer multiple of period of metasurface)
    k=1.1k₀ # wavenumber of incident field
    n = 10 # number of modes to scatter
    dims = ((n, L), (1, 1.0))
    I₁ = [i == j == 1 ? 1.0 : 0.0 for i in 1:n, j in 1:2]
    I₂ = zeros(n, 2)
    gradient(s -> cost(solve(DeltaRCWAProblem(s, dims, k, I₁, I₂))), sheet)
end # ((θ = -0.0, k = 0.0, d = 0.0),)