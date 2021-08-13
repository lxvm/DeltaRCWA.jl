export UniformSlab

"""
    UniformSlab(depth::T where T <: Real)

Wrapper for a real number needed to propagate PlanewaveModes that distance
"""
struct UniformSlab{T} <: RCWASlab{T, 0}
    depth::T
end

"""
    smatrix(slab::UniformSlab, modes::PlanewaveModes)

Returns the scattering matrix due to mode propagation through the z-direction by
adding a phase to each mode corresponding to the wavevector component kz, its
sign, and the depth of the layer
"""
function smatrix(slab::UniformSlab, modes, ::AbstractPolarization)::BlockMatrix
    # add phase from propagating the plane wave in the z direction through layer
    # fw means forward: from lower to larger z
	fwpropagator = Matrix(Diagonal(exp.(im * modes.kz * slab.depth)))
    # bk means backward: from larger to lower z
    bkpropagator = conj.(fwpropagator)
    mortar(
        (zeros(size(fwpropagator)), bkpropagator),
        (fwpropagator,              zeros(size(bkpropagator))),
    )
end