# Defines the trivial fall-back methods that should be implemented by `RCWASheet`s
export σₑˣˣ, σₑˣʸ, σₑʸˣ, σₑʸʸ, σₘˣˣ, σₘˣʸ, σₘʸˣ, σₘʸʸ

nonconducting(x⃗) = zeros(Bool, length.(x⃗))

σₑˣˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₑˣʸ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₑʸˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₑʸʸ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘˣˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘˣʸ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘʸˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘʸʸ(::RCWASheet, x⃗) = nonconducting(x⃗)

"""
    smatrix(sheet::RCWASheet{1}, modes, ::UncoupledPolarization)

Returns a 2x2 BlockMatrix for the scattering of modes specific to the TE or TM 
polarization
"""
function smatrix(sheet::RCWASheet{T, 1} where T, modes, pol::UncoupledPolarization)
    n = length(modes.kz)
    BlockMatrix(
        _1Dsheetsmatrix(_get_params_1Dsheetsmatrix(sheet, modes, pol)...),
        [n, n], [n, n]
    )
end

_get_params_1Dsheetsmatrix(sheet, modes, ::TE) = (
    σₘˣˣ(sheet, modes.x⃗), σₑʸʸ(sheet, modes.x⃗), modes.kz, modes.ω * modes.M.μ,
)
_get_params_1Dsheetsmatrix(sheet, modes, ::TM) = (
    σₑˣˣ(sheet, modes.x⃗), σₘʸʸ(sheet, modes.x⃗), modes.kz, modes.ω * modes.M.ϵ,
)

"construct the dense scattering matrix"
function _1Dsheetsmatrix(σˣˣ, σʸʸ, kz, pol_ω)
    kz_term = Diagonal(reshape(kz, :)) / pol_ω
    σ̃ˣˣ_term, σ̃ʸʸ_term = diagFT.((0.5σˣˣ, 0.5σʸʸ))
    A = _get_matrix_scattered_side(σ̃ˣˣ_term, σ̃ʸʸ_term, kz_term)
    B = _get_matrix_incident_side(σ̃ˣˣ_term, σ̃ʸʸ_term, kz_term)
    A\B
end

"""
    diagFT(A::AbstractArray{<:Number, N}) where N

Accepts a diagonal endomorphism in N-D position space and returns the Fourier
transformed  rank-2N tensor reshaped to a matrix with flattened indices
"""
function diagFT(A::AbstractArray{<:Number, ndim}) where ndim
    n⃗ = size(A)
    n = prod(n⃗)
    Ã = Matrix(reshape(Diagonal(reshape(A, n)), (n⃗..., n⃗...)))
    Ã = fft(ifft(Ã, 1:ndim), (ndim+1):2ndim)
    reshape(Ã, (n, n))
end

function _get_matrix_scattered_side(σ̃ˣˣ_term, σ̃ʸʸ_term, kz_term)
    [
        (I + σ̃ˣˣ_term * kz_term)      (-I - σ̃ˣˣ_term * kz_term);
        (kz_term + σ̃ʸʸ_term)                (kz_term + σ̃ʸʸ_term)
    ]
end

function _get_matrix_incident_side(σ̃ˣˣ_term, σ̃ʸʸ_term, kz_term)
    [
        (-I + σ̃ˣˣ_term * kz_term)    (I - σ̃ˣˣ_term * kz_term);
        (kz_term - σ̃ʸʸ_term)           (kz_term - σ̃ʸʸ_term)
    ]
end

"""
    smatrixLinearMap(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a RCWASheet in a matrix-free fashion.
"""
function smatrixLinearMap(sheet::RCWASheet{T, 1} where T, modes, pol::UncoupledPolarization)
    _1DsheetsmatrixLinearMap(_get_params_1Dsheetsmatrix(sheet, modes, pol)...)
end

function _1DsheetsmatrixLinearMap(σˣˣ, σʸʸ, kz, pol_ω)
    kz_term = Diagonal(reshape(kz, :)) / pol_ω
    σˣˣ_term, σʸʸ_term = Diagonal.(reshape.((0.5σˣˣ, 0.5σʸʸ), :))
    A = _get_LinearMap_scattered_side(σˣˣ_term, σʸʸ_term, kz_term)
    B = _get_LinearMap_incident_side(σˣˣ_term, σʸʸ_term, kz_term)
    LinearMap(x -> gmres(A, B(x)), 2size(kz_term, 1))
end

function _get_LinearMap_scattered_side(σˣˣ_term, σʸʸ_term, kz_term)
    N = size(kz_term, 1)
    LinearMap(
        function (x::AbstractVector)
            I₁ = x[1:N]
            I₂ = x[(N+1):2N]
            sumI₁I₂ = I₁ + I₂
            diffI₁I₂ = I₁ - I₂
            vcat(
                diffI₁I₂ + ifft(σˣˣ_term * fft(kz_term * diffI₁I₂)),
                (kz_term * sumI₁I₂) + ifft(σʸʸ_term * fft(sumI₁I₂))
            )
        end,
        2N
    )
end

function _get_LinearMap_incident_side(σˣˣ_term, σʸʸ_term, kz_term)
    N = size(kz_term, 1)
    LinearMap(
        function (x::AbstractVector)
            I₁ = x[1:N]
            I₂ = x[(N+1):2N]
            sumI₁I₂ = I₁ + I₂
            diffI₁I₂ = I₁ - I₂
            vcat(
                -diffI₁I₂ + ifft(σˣˣ_term * fft(kz_term * diffI₁I₂)),
                (kz_term * sumI₁I₂) - ifft(σʸʸ_term * fft(sumI₁I₂))
            )
        end,
        2N
    )
end

"""
    smatrixBlockMap(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a RCWASheet in a matrix-free fashion.
"""
function smatrixBlockMap(sheet::RCWASheet{T, 1} where T, modes, pol::UncoupledPolarization)
    _1DsheetsmatrixBlockMap(_get_params_1Dsheetsmatrix(sheet, modes, pol)...)
end

function _1DsheetsmatrixBlockMap(σˣˣ, σʸʸ, kz, pol_ω)::BlockMap
    N = length(kz)
    kz_term = LinearMap(Diagonal(reshape(kz, :)) / pol_ω)
    σˣˣ_term, σʸʸ_term = LinearMap.(Diagonal.(reshape.((0.5σˣˣ, 0.5σʸʸ), :)))
    _fft = LinearMap(x -> fft(x), N)
    _ifft = LinearMap(x -> ifft(x), N)
    A = _get_BlockMap_scattered_side(_fft, _ifft, σˣˣ_term, σʸʸ_term, kz_term)
    # A = _get_LinearMap_scattered_side(σˣˣ_term, σʸʸ_term, kz_term)
    B = _get_BlockMap_incident_side(_fft, _ifft, σˣˣ_term, σʸʸ_term, kz_term)
    invA₁₁ = LinearMap(x -> gmres(A, vcat(x, zero(x)))[1:N], N)
    invA₁₂ = LinearMap(x -> gmres(A, vcat(zero(x), x))[1:N], N)
    invA₂₁ = LinearMap(x -> gmres(A, vcat(x, zero(x)))[(N+1):2N], N)
    invA₂₂ = LinearMap(x -> gmres(A, vcat(zero(x), x))[(N+1):2N], N)
    B₁₁, B₁₂, B₂₁, B₂₂ = B.maps
    [
        invA₁₁ * B₁₁ + invA₁₂ * B₂₁     invA₁₁ * B₁₂ + invA₁₂ * B₂₂;
        invA₂₁ * B₁₁ + invA₂₂ * B₂₁     invA₂₁ * B₁₂ + invA₂₂ * B₂₂
    ]
end

function _get_BlockMap_scattered_side(_fft, _ifft, σˣˣ_term, σʸʸ_term, kz_term)
    [
        I + _ifft * σˣˣ_term * _fft * kz_term     -I - _ifft * σˣˣ_term * _fft * kz_term;
        kz_term + _ifft * σʸʸ_term * _fft         kz_term + _ifft * σʸʸ_term * _fft
    ]
end

function _get_BlockMap_incident_side(_fft, _ifft, σˣˣ_term, σʸʸ_term, kz_term)
    [
        -I + _ifft * σˣˣ_term * _fft * kz_term    I + _ifft * σˣˣ_term * _fft * kz_term;
        kz_term - _ifft * σʸʸ_term * _fft         kz_term - _ifft * σʸʸ_term * _fft
    ]
end

"""
    smatrix(sheet::RCWASheet{2}, modes, ::UncoupledPolarization)

Returns a 2x2 BlockMatrix for the scattering of modes specific to the TE or TM 
polarization
"""
function smatrix(sheet::RCWASheet{T, 2} where T, modes, ::CoupledPolarization)
    n = length(modes.kz)
    BlockMatrix(
        _2Dsheetsmatrix(_get_params_2Dsheetsmatrix(sheet, modes)...),
        [2n, 2n], [2n, 2n]
    )
end

"Return the variables needed to construct the dense scattering matrix"
function _get_params_2Dsheetsmatrix(sheet, modes)
    n = length(modes.kz)
    R = (
        xx = (0I)(n),
        xy = (-I)(n),
        yx = (1I)(n),
        yy = (0I)(n),
    )
    ωμ = modes.ω * modes.M.μ
    k⃗² = modes.ω * modes.M.ϵ * ωμ
    ρₑ = (
        xx = ρₑˣˣ(sheet, modes.x⃗),
        xy = ρₑˣʸ(sheet, modes.x⃗),
        yx = ρₑʸˣ(sheet, modes.x⃗),
        yy = ρₑʸʸ(sheet, modes.x⃗),
    )
    σₘ = (
        xx = σₘˣˣ(sheet, modes.x⃗),
        xy = σₘˣʸ(sheet, modes.x⃗),
        yx = σₘʸˣ(sheet, modes.x⃗),
        yy = σₘʸʸ(sheet, modes.x⃗),
    )
    K = (
        xx = [k⃗[1]^2 for k⃗ in Iterators.product(modes.k⃗)],
        xy = [k⃗[1] * k⃗[2] for k⃗ in Iterators.product(modes.k⃗)],
        # yx = [k⃗[2] * k⃗[1] for k⃗ in Iterators.product(modes.k⃗)],
        yy = [k⃗[2]^2 for k⃗ in Iterators.product(modes.k⃗)],
    )
    ρₑ, σₘ, R, K, modes.kz, k⃗², ωμ
end

function _2Dsheetsmatrix(ρₑ, σₘ, R, K, kz, k⃗², ωμ)::BlockMatrix
    K_term = mortar(map(x -> @. Diagonal(reshape(x, :)), [
        (k⃗² .+ K.xx, K.xy),
        (K.xy, k⃗² .+ K.yy),
    ])...)
    R_term = mortar(
        (R.xx, R.xy),
        (R.yx, R.yy),
    )
    # need to do Fourier transform of each block
    ρ̃_term = mortar(map(x -> @. diagFT(x), [
        (ρₑ.xx, ρₑ.xy),
        (ρₑ.yx, ρₑ.yy),
    ])...)
    ρ̃_term = Matrix(K_term * transpose(R_term) * ρ̃_term * R_term)
    σ̃_term = mortar(map(x -> @. diagFT(x), [
        (σₘ.xx, σₘ.xy),
        (σₘ.yx, σₘ.yy),
    ])...)
    σ̃_term = Matrix(K_term * σ̃_term)
        kz_term = ωμ * Diagonal(vcat(reshape.((kz, kz), :)...))
    A = _get_2Dsmatrix_scattered_side(ρ̃_term, σ̃_term, kz_term)
    B = _get_2Dsmatrix_incident_side(ρ̃_term, σ̃_term, kz_term)
    A\B
end

function _get_2Dsmatrix_scattered_side(ρ̃_term, σ̃_term, kz_term)
    [
        ρ̃_term - 0.5*kz_term   -ρ̃_term + 0.5*kz_term;
        σ̃_term + 2*kz_term     σ̃_term + 2*kz_term
    ]
end

function _get_2Dsmatrix_incident_side(ρ̃_term, σ̃_term, kz_term)
    [
        ρ̃_term + 0.5*kz_term   -ρ̃_term - 0.5*kz_term;
        σ̃_term - 2*kz_term     σ̃_term - 2*kz_term
    ]
end

function smatrixLinearMap end
function _2DsheetsmatrixLinearMap end

function smatrixBlockMap end
function _2DsheetsmatrixBlockMap end