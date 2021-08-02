export UniformSlab

"""
    UniformSlab(depth::Real, ϵ::Number, μ::Number)

Stores data needed to propagate a wave in a uniform medium
"""
struct UniformSlab <: RCWASlab{0}
    depth::T where T <: Real
    M::UniformMedium
end

"""
    smatrix(slab::UniformSlab, modes::PlanewaveModes)

Returns the scattering matrix due to mode propagation through the z-direction by
adding a phase to each mode corresponding to the wavevector component kz, its
sign, and the depth of the layer
"""
function smatrix(slab::UniformSlab, modes::PlanewaveModes)::BlockMatrix
    # @assert slab.M == modes.M
    # add phase from propagating the plane wave in the z direction through layer
    # fw means forward: from lower to larger z
	fwpropagator = Matrix(Diagonal(exp.(im * modes.kz * slab.depth)))
    # bk means backward: from larger to lower z
    bkpropagator = conj.(fwpropagator)
    return mortar(reshape([
        zeros(size(fwpropagator)),
		fwpropagator,
        bkpropagator,
		zeros(size(bkpropagator)),
    ], 2, 2))
end