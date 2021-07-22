"""
Stores data needed to model a scattering interface of 2 uniform media
"""
struct UniformInterface <: DeltaScatterer
    ϵ₁::Number
    μ₁::Number
    ϵ₂::Number
    μ₂::Number
end