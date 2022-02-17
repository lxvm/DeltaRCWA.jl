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

const θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction
const θᵗ = π/3 # desired transmitted field angle
const d = cos(θᵗ)-cos(θⁱ) # design parameter
const k₀ = 10.0 # wavenumber that attains the metasurface design goals
const λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals
const L₀ = λ₀/abs(d) # period of metasurface
const sheet = ComplexExpSheet(θⁱ, k₀, d) # sheet to scatter

const L=L₀ # length of domain (some integer multiple of period of metasurface)
const k=1.1k₀ # wavenumber of incident field
const n = 20 # number of modes to scatter
const dims = ((n, L),)
const I₁ = [i == 1 ? 1.0 : 0.0 for i in 1:n]
const I₂ = zeros(n)


function cost(sol::DeltaRCWA.DeltaRCWASolution{1})
    norm(sol.O₂) + norm(sol.O₁)
end

ans = gradient(s -> cost(solve(DeltaRCWAProblem(s, dims, k, I₁, I₂))), sheet)

# MWE of error: cannot broadcast over a product Iterator
# julia> gradient(k -> norm(sum.(Iterators.product(k...))), (1:10,))