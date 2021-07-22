"""
    UniformLayer(depth::Real, ϵ::Number, μ::Number)

Stores data needed to propagate a wave in a uniform medium
"""
struct UniformLayer <: AbstractScatterer
    depth::Real
    ϵ::Number
    μ::Number
end

Air(depth) = UniformLayer(depth, 1, 1)

"""
    get_kz(layer::UniformLayer, k::NTuple{N, Frequencies} where N, ω::Real)

Returns a vector of kz from the dispersion relation k⋅k = ω²ϵμ
"""
function get_kz(layer::UniformLayer, k::Tuple{}, ω::Real)
	return [sqrt(Complex((ω^2) * layer.ϵ * layer.μ))]
end

function get_kz(layer::UniformLayer, k::NTuple{N, Frequencies} where N, ω::Real)
    λ⁻² = (ω^2) * layer.ϵ * layer.μ
	return [sqrt(Complex(λ - sum(abs2, e))) for e in Iterators.product(k...)]
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
	fwpropagator = Matrix(Diagonal(exp.(im* get_kz(layer, k, ω) * layer.depth)))
    # bk means backward: from larger to lower z
    bkpropagator = conj.(fwpropagator)
    return mortar(reshape([
        zeros(size(fwpropagator)),
		fwpropagator,
        bkpropagator,
		zeros(size(bkpropagator)),
    ], 2, 2))
end