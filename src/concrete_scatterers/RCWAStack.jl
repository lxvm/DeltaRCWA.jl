export RCWAStack, UniformMediumSheetStack

struct RCWAStack{N} <: RCWAScatterer{N, 3}
    layers::NTuple{M, RCWAScatterer} where M
end

"""
    RCWAStack(layers::AbstractVector{AbstractScatterer})

Creates a scatterer that is a sequence of other scatterers
TODO: this stack needs to carry information about the medium! It should be a
1D graph whose vertices are of type RCWASlab and edges are of type RCWASheet
"""
function RCWAStack(layers::T) where T <: AbstractVector{RCWAScatterer}
    N, M = find_largest_dim(layers)
    @assert M == 3
    RCWAStack{N}(layers)
end

"Return the largest value of N for which any element is <: RCWAScatterer{N}"
function find_largest_dim(layers::T) where T <: AbstractVector{RCWAScatterer}
    M = maxMdim
    while M >= minMdim
        N = maxNdim
        while N >= minNdim
            for layer in layers
                if typeof(layer) <: RCWAScatterer{N, M}
                    return N, M
                end
            end
            N -= 1
        end
        M -= 1
    end
    return N, M
end

"""
    smatrix(stack::RCWAStack, modes::PlanewaveModes)

Computes the full scattering matrix of the stack of scatterers using the
[Redheffer star product](https://en.wikipedia.org/wiki/Redheffer_star_product).
"""
function smatrix(stack::RCWAStack, modes::PlanewaveModes)
    S = smatrix(stack.layers[1], modes)
    for layer in stack.layers[2:end]
        S = smat_star(S, smatrix(layer, modes))
    end
    return S
end

"""
A simple tuple of sheets with a complementary tuple of the spacings between
sheets in a UniformMedium
"""
struct UniformMediumSheetStack{N, L, D, T} <: RCWAScatterer{N, 3}
    sheets::NTuple{L, RCWASheet{N}} # this is a tuple of pointers
    depths::NTuple{D, <: Real} # this is a tuple of a concrete subtype of real
    medium::UniformMedium{T}
    function UniformMediumSheetStack(
        sheets::NTuple{L, RCWASheet{N}},
        depths::NTuple{D, <: Real},
        medium::UniformMedium{T},
    ) where {N, L, D, T}
        L-1==D || throw(ArgumentError("there must be one gap less than the number of sheets"))
        new{N, L, D, T}(sheets, depths, medium)
    end
end