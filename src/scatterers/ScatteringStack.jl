"""
    ScatteringStack(layers::AbstractVector{AbstractScatterer})

Creates a scatterer that is a sequence of other scatterers
"""
mutable struct ScatteringStack <: AbstractScatterer
    layers::AbstractVector{AbstractScatterer}
end

"""
    smatrix(stack::ScatteringStack, k::NTuple{N, Frequencies} where N, ω::Real)

Computes the full scattering matrix of the stack of scatterers using the
[Redheffer star product](https://en.wikipedia.org/wiki/Redheffer_star_product).
"""
function smatrix(stack::ScatteringStack, k::NTuple{N, Frequencies} where N, ω::Real)
    S = smatrix(stack.layers[1], k, ω)
    for layer in stack.layers[2:end]
        S = smat_star(S, smatrix(layer, k, ω))
    end
    return S
end