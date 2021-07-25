"""
Stores data needed to model a scattering interface of 2 uniform media
"""
struct UniformInterface{T <: AbstractUniformMedium} <: DeltaScatterer
    medium₁::T
    medium₂::T
end

# passing through an interface like this one requires that the kx, ky
# components of the wavevector remain the same, but the new wave velocity
# affects the value of kz in the medium, so the mapping between PlanewaveModes
# is one-to-one on the product space of (kx, ky)