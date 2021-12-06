export UniformSlab

"""
    UniformSlab(depth::Float64)

Wrapper for a `Float64` needed to propagate PlaneWaves that distance
"""
struct UniformSlab
    depth::Float64
end

"""
    _get_propagators(kz::AbstractArray, Δz::Real; reps=1)

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

function _sMatrix(pw::PlaneWaves{N}, slab::UniformSlab, um::UniformMedium) where N
    kz = _get_kz(pw, um)
    fw, bk = _get_propagators(kz, slab.depth; reps=N)
    zo = (Bool(0) * I)(N*length(kz))
    [
        zo bk;
        fw zo;
    ]
end

function _sBlockMatrix(pw::PlaneWaves{N}, slab::UniformSlab, um::UniformMedium)::BlockMatrix where N
    kz = _get_kz(pw, um)
    fw, bk = _get_propagators(kz, slab.depth; reps=N)
    zo = (Bool(0) * I)(N*length(kz))
    mortar(
        (zo, bk),
        (fw, zo),
    )
end

function _sBlockMap(pw::PlaneWaves{N}, slab::UniformSlab, um::UniformMedium)::BlockMap where N
    kz = _get_kz(pw, um)
    fw, bk = _get_propagators(kz, slab.depth; reps=N)
    [
        Bool(0) * I    LinearMap(bk);
        LinearMap(fw)    Bool(0) * I
    ]
end

function _sLinearMap(pw::PlaneWaves{N}, slab::UniformSlab, um::UniformMedium)::LinearMap where N
    kz = _get_kz(pw, um)
    d = length(kz)
    fw, bk = _get_propagators(kz, slab.depth; reps=N)
    LinearMap(2d*N) do x
        vcat(
            bk * x[(d*N+1):2d*N],
            fw * x[1:d*N],
        )
    end
end