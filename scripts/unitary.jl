using DeltaRCWA
using LinearAlgebra
using Statistics

struct ComplexExpSheet{T} <: Sheet
    θ::T # incidence angle
    k::T # wavenumber
    d::T # designer parameter
end

function DeltaRCWA.Zₑˣˣ(sheet::ComplexExpSheet, x⃗::Tuple)
    z = -0.5sin(sheet.θ)*(1-exp(1im*sheet.k*sheet.d*x⃗[1]))
    im*imag(z)
end

function DeltaRCWA.Yₘʸʸ(sheet::ComplexExpSheet, x⃗::Tuple)
    y = -2sin(sheet.θ)*(1+exp(1im*sheet.k*sheet.d*x⃗[1]))
    im*imag(y)
end