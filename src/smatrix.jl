# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers

function smatrix(layer <: AbstractLayer, wavenumbers <: AbstractVector)
    # currently assuming normal incidence, so the wavenumbers are due to the
    # kz component
    error("NotImplemented")
end

function smatrix(layer <: UniformLayer, wavenumbers <: AbstractVector, ω₀)
    propagator = @. Diagonal(Complex(exp(
        im * sqrt(Complex(ω₀ * layer.ϵ * layer.μ - wavenumbers^2)) * layer.depth
    )))
    return mortar([
        zeros(size(propagator)) propagator;
        propagator zeros(size(propagator))
    ])
end
