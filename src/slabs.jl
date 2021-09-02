export UniformSlab

"""
    UniformSlab(depth::T)

Wrapper for a real number needed to propagate PlanewaveModes that distance
"""
struct UniformSlab{T<:Real} <: RCWASlab{0}
    depth::T
end

"""
    _get_propagators(kz::AbstractVector, Δz::Real)

Calculate phases to propagate the plane wave in the +/- z direction
"""
function _get_propagators(kz::AbstractVector, Δz::Real)
    # fw means forward: from lower to larger z (port i to i+1)
    fw = Diagonal(exp.((im * Δz) * kz))
    # bk means backward: from larger to lower z (port i+1 to i)
    bk = conj.(fw)
    fw, bk
end

"""
    smatrix(slab::UniformSlab, modes::PlanewaveModes)

Returns the scattering matrix due to mode propagation through the z-direction by
adding a phase to each mode corresponding to the wavevector component kz, its
sign, and the depth of the layer
"""
function smatrix(slab::UniformSlab, modes, pol)::BlockMatrix
    fw, bk = _get_propagators(reshape(modes.kz, :), slab.depth)
    zo = (Bool(0) * I)(length(modes.kz))
    mortar(
        (zo, bk),
        (fw, zo),
    )
end

function smatrixBlockMap(slab::UniformSlab, modes, pol)::BlockMap
    fw, bk = _get_propagators(reshape(modes.kz, :), slab.depth)
    [
        Bool(0) * I    LinearMap(bk);
        LinearMap(fw)    Bool(0) * I
    ]
end

function smatrixLinearMap(slab::UniformSlab, modes, pol)::LinearMap
    d = length(modes.kz)
    fw, bk = _get_propagators(reshape(modes.kz, :), slab.depth)
    LinearMap(x -> vcat(
            bk .* x[(d+1):2d],
            fw .* x[1:d],
        ),
        2d
    )
end