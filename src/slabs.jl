export UniformSlab

"""
    UniformSlab(depth::T)

Wrapper for a real number needed to propagate PlanewaveModes that distance
"""
struct UniformSlab{T<:Real}
    depth::T
end

"""
    _get_propagators(kz::AbstractVector, Δz::Real; reps=1)

Calculate phases to propagate the plane wave in the +/- z direction.
Keyword arguments `reps` repeats the vector several times
"""
function _get_propagators(kz::AbstractArray, Δz::Real; reps=1)
    # fw means forward: from lower to larger z (port i to i+1)
    fw = repeat(exp.((im * Δz) * vec(kz)), reps)
    # bk means backward: from larger to lower z (port i+1 to i)
    bk = conj.(fw)
    Diagonal.((fw, bk))
end

"""
    smatrix(slab::UniformSlab, modes::PlanewaveModes)

Returns the scattering matrix due to mode propagation through the z-direction by
adding a phase to each mode corresponding to the wavevector component kz, its
sign, and the depth of the layer
"""
function _sMatrix(slab::UniformSlab, modes::PlanewaveModes{T, N}, pol) where {T, N}
    fw, bk = _get_propagators(modes.kz, slab.depth; reps=N)
    zo = (Bool(0) * I)(N*length(modes.kz))
    [
        zo bk;
        fw zo;
    ]
end

function _sBlockMatrix(slab::UniformSlab, modes::PlanewaveModes{T, N}, pol)::BlockMatrix where {T, N}
    fw, bk = _get_propagators(modes.kz, slab.depth; reps=N)
    zo = (Bool(0) * I)(N*length(modes.kz))
    mortar(
        (zo, bk),
        (fw, zo),
    )
end

function _sBlockMap(slab::UniformSlab, modes::PlanewaveModes{T, N}, pol)::BlockMap where {T, N}
    fw, bk = _get_propagators(modes.kz, slab.depth; reps=N)
    [
        Bool(0) * I    LinearMap(bk);
        LinearMap(fw)    Bool(0) * I
    ]
end

function _sLinearMap(slab::UniformSlab, modes::PlanewaveModes{T, N}, pol)::LinearMap where {T, N}
    d = length(modes.kz)
    fw, bk = _get_propagators(modes.kz, slab.depth; reps=N)
    LinearMap(x -> vcat(
            bk * x[(d*N+1):2d*N],
            fw * x[1:d*N],
        ),
        2d*N
    )
end