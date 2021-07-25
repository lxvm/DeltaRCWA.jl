"""
    ScatteringStack(layers::AbstractVector{AbstractScatterer})

Creates a scatterer that is a sequence of other scatterers
"""
mutable struct ScatteringStack{T <: AbstractVector} <: AbstractScatterer
    layers::T
end

"""
    smatrix(stack::ScatteringStack, k::NTuple{N, Frequencies} where N, ω::Real)

Computes the full scattering matrix of the stack of scatterers using the
[Redheffer star product](https://en.wikipedia.org/wiki/Redheffer_star_product).
"""
function smatrix(stack::ScatteringStack, modes)
    S = smatrix(stack.layers[1], modes)
    for layer in stack.layers[2:end]
        S = smat_star(S, smatrix(layer, modes))
    end
    return S
end