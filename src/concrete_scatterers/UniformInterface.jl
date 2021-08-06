export UniformInterface

"""
Stores data needed to model a scattering interface of 2 uniform media
"""
struct UniformInterface <: RCWASheet{0}
    M₁::UniformMedium
    M₂::UniformMedium
end

# passing through an interface like this one requires that the kx, ky
# components of the wavevector remain the same, but the new wave velocity
# affects the value of kz in the medium, so the mapping between PlanewaveModes
# is one-to-one on the product space of (kx, ky)

"""
    smatrix(inf::UniformInterface, k₁::PlanewaveModes, k₂::PlanewaveModes)

Uses Fresnel coefficients to calculate the transmittance and reflectance of a
periodic set of modes and forms them into a scattering matrix.
"""
function smatrix(inf::UniformInterface, k₁::PlanewaveModes, k₂::PlanewaveModes)
    # port 1 -> 2
    β₂₁ = sqrt((inf.M₁.μ * inf.M₂.ϵ) / (inf.M₁.ϵ * inf.M₂.μ))
    n₂₁ = sqrt((inf.M₁.ϵ * inf.M₁.μ) / (inf.M₂.ϵ * inf.M₂.μ))
    α₂₁ = @. n₂₁ * sqrt(k₂.kz / k₁.kz)
    # port 2 -> 1
    β₁₂ = 1 / β₂₁
    α₁₂ = @. 1 / α₂₁
    # Fresnel coefficients by mode
    S₁₁ = Diagonal(@. (α₂₁ - β₂₁) / (α₂₁ + β₂₁)) # R₁
    S₂₁ = Diagonal(@. 2 / (α₂₁ + β₂₁)) # T₁
    S₁₂ = Diagonal(@. 2 / (α₁₂ + β₁₂)) # T₂
    S₂₂ = Diagonal(@. (α₁₂ - β₁₂) / (α₁₂ + β₁₂)) # R₂
    mortar(reshape([S₁₁, S₂₁, S₁₂, S₂₂], 2, 2))
end
