"""
    UniformLayer(depth::Real, ϵ::Number, μ::Number)

Stores data needed to propagate a wave in a uniform medium
"""
struct UniformLayer <: DeltaScatterer
    depth::Real
    medium::AbstractUniformMedium
end

"""
    smatrix(layer::UniformLayer, k::NTuple{N, Frequencies} where N, ω::Real)

Returns the scattering matrix due to mode propagation through the z-direction by
adding a phase to each mode corresponding to the wavevector component kz, its
sign, and the depth of the layer
"""
function smatrix(layer::UniformLayer, k::NTuple{N, Frequencies} where N, ω::Real)
    # add phase from propagating the plane wave in the z direction through layer
    # fw means forward: from lower to larger z
	fwpropagator = Matrix(Diagonal(exp.(im * get_kz(layer.medium, k, ω) * layer.depth)))
    # bk means backward: from larger to lower z
    bkpropagator = conj.(fwpropagator)
    return mortar(reshape([
        zeros(size(fwpropagator)),
		fwpropagator,
        bkpropagator,
		zeros(size(bkpropagator)),
    ], 2, 2))
end