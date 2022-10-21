using DeltaRCWA
using LinearAlgebra

function _get_wave_params(ω, kx, ϵ, μ)
    kz = sqrt(Complex(ω^2*ϵ*μ - kx^2))
    Y = sqrt(Complex(ϵ/μ))
    v = sqrt(Complex(ϵ*μ))^-1
    kz, Y, v
end

function FresnelTM(ω, kx, ϵ₁, μ₁, ϵ₂, μ₂)
    kz₁, Y₁, v₁ = _get_wave_params(ω, kx, ϵ₁, μ₁)
    kz₂, Y₂, v₂ = _get_wave_params(ω, kx, ϵ₂, μ₂)
    S = [
        -(Y₁*v₂*kz₂ - Y₂*v₁*kz₁)   2*Y₁*v₂*kz₂
        2*Y₂*v₁*kz₁   -(Y₂*v₁*kz₁ - Y₁*v₂*kz₂)
        # Y₁*v₂*kz₂ - Y₂*v₁*kz₁   2*Y₂*v₂*kz₂
        # 2*Y₁*v₁*kz₁   Y₂*v₁*kz₁ - Y₁*v₂*kz₂
    ]
    return S / (Y₁*v₂*kz₂ + Y₂*v₁*kz₁)
end
function uFresnelTM(ω, kx, ϵ₁, μ₁, ϵ₂, μ₂)
    kz₁, Y₁, v₁ = _get_wave_params(ω, kx, ϵ₁, μ₁)
    kz₂, Y₂, v₂ = _get_wave_params(ω, kx, ϵ₂, μ₂)
    S = FresnelTM(ω, kx, ϵ₁, μ₁, ϵ₂, μ₂)
    α = sqrt(sqrt((Y₂*v₁*kz₁) / (Y₁*v₂*kz₂)))
    A = Diagonal([α^-1, α])
    inv(A)*S*A
end
#=
function FresnelTE(ω, kx, ϵ₁, μ₁, ϵ₂, μ₂)
    kz₁, Y₁, v₁ = _get_wave_params(ω, kx, ϵ₁, μ₁)
    kz₂, Y₂, v₂ = _get_wave_params(ω, kx, ϵ₂, μ₂)
    S = [
        Y₁*v₁*kz₁ - Y₂*v₂*kz₂   2*Y₂*v₂*kz₂
        2*Y₁*v₁*kz₁   Y₂*v₂*kz₂ - Y₁*v₁*kz₁
    ]
    return S / (Y₁*v₁*kz₁ + Y₂*v₂*kz₂)
end
function uFresnelTE(ω, kx, ϵ₁, μ₁, ϵ₂, μ₂)
    kz₁, Y₁, v₁ = _get_wave_params(ω, kx, ϵ₁, μ₁)
    kz₂, Y₂, v₂ = _get_wave_params(ω, kx, ϵ₂, μ₂)
    S = FresnelTE(ω, kx, ϵ₁, μ₁, ϵ₂, μ₂)
    α = sqrt(sqrt((Y₁*v₁*kz₁) / (Y₂*v₂*kz₂)))
    A = Diagonal([α^-1, α])
    inv(A)*S*A
end
=#
struct TransparentSheet <: Sheet end
DeltaRCWA.MagneticResponseStyle(::Type{TransparentSheet}) = Admittance()
DeltaRCWA.ElectricResponseStyle(::Type{TransparentSheet}) = Admittance()

function blockify(i, n)
    d, r = divrem(i-1, n)
    2*r+1+d # in general replace 2 with n_block
end