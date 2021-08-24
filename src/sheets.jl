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
        _1Dsheetsmatrix(_get_params_smatrix(sheet, modes, pol)...),
        [n, n], [n, n]
    )
end

"convenience function to solve problems for each polarization"
_getσωbypol(sheet, modes, ::TE) = (modes.ω * modes.M.μ, σₘˣˣ(sheet, modes.x⃗), σₑʸʸ(sheet, modes.x⃗))
_getσωbypol(sheet, modes, ::TM) = (modes.ω * modes.M.ϵ, σₑˣˣ(sheet, modes.x⃗), σₘʸʸ(sheet, modes.x⃗))

function _get_params_smatrix(sheet, modes, pol)
    pol_ω, σˣˣ, σʸʸ = _getσωbypol(sheet, modes, pol)
    halfσˣˣ = 0.5Diagonal(reshape(σˣˣ, :))
    halfσʸʸ = 0.5Diagonal(reshape(σʸʸ, :))
    pol_kz = Diagonal(reshape(modes.kz, :)) / pol_ω
    halfσˣˣ, halfσʸʸ, pol_kz
end

"construct the dense scattering matrix"
function _1Dsheetsmatrix(halfσˣˣ::Diagonal, halfσʸʸ::Diagonal, pol_kz::Diagonal)
    # tilde implies Fourier basis
    halfσ̃ˣˣ = ifft(fft(Matrix(halfσˣˣ), 2), 1)
    halfσ̃ʸʸ = ifft(fft(Matrix(halfσʸʸ), 2), 1)
    A = _get_matrix_scattered_side(halfσ̃ˣˣ, halfσ̃ʸʸ, pol_kz)
    B = _get_matrix_incident_side(halfσ̃ˣˣ, halfσ̃ʸʸ, pol_kz)
    A\B
end

function _get_matrix_scattered_side(halfσ̃ˣˣ, halfσ̃ʸʸ, pol_kz)
    [
        (I + halfσ̃ˣˣ * pol_kz)      (-I - halfσ̃ˣˣ * pol_kz);
        (pol_kz + halfσ̃ʸʸ)                (pol_kz + halfσ̃ʸʸ)
    ]
end

function _get_matrix_incident_side(halfσ̃ˣˣ, halfσ̃ʸʸ, pol_kz)
    [
        (-I + halfσ̃ˣˣ * pol_kz)    (I - halfσ̃ˣˣ * pol_kz);
        (pol_kz - halfσ̃ʸʸ)           (pol_kz - halfσ̃ʸʸ)
    ]
end

"""
    smatrixLinearMap(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a RCWASheet in a matrix-free fashion.
"""
function smatrixLinearMap(sheet::RCWASheet{T, 1} where T, modes, pol::UncoupledPolarization)
    _1DsheetsmatrixLinearMap(_get_params_smatrix(sheet, modes, pol)...)
end

function _1DsheetsmatrixLinearMap(halfσˣˣ::Diagonal, halfσʸʸ::Diagonal, pol_kz::Diagonal)
    A = _get_LinearMap_scattered_side(halfσˣˣ, halfσʸʸ, pol_kz)
    B = _get_LinearMap_incident_side(halfσˣˣ, halfσʸʸ, pol_kz)
    LinearMap(x -> gmres(A, B(x)), 2size(pol_kz, 1))
end

function _get_LinearMap_scattered_side(halfσˣˣ, halfσʸʸ, pol_kz)
    N = size(pol_kz, 1)
    LinearMap(
        function (x::AbstractVector)
            I₁ = x[1:N]
            I₂ = x[(N+1):2N]
            sumI₁I₂ = I₁ + I₂
            diffI₁I₂ = I₁ - I₂
            vcat(
                diffI₁I₂ + ifft(halfσˣˣ * fft(pol_kz * diffI₁I₂)),
                (pol_kz * sumI₁I₂) + ifft(halfσʸʸ * fft(sumI₁I₂))
            )
        end,
        2N
    )
end

function _get_LinearMap_incident_side(halfσˣˣ, halfσʸʸ, pol_kz)
    N = size(pol_kz, 1)
    LinearMap(
        function (x::AbstractVector)
            I₁ = x[1:N]
            I₂ = x[(N+1):2N]
            sumI₁I₂ = I₁ + I₂
            diffI₁I₂ = I₁ - I₂
            vcat(
                -diffI₁I₂ + ifft(halfσˣˣ * fft(pol_kz * diffI₁I₂)),
                (pol_kz * sumI₁I₂) - ifft(halfσʸʸ * fft(sumI₁I₂))
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
    _1DsheetsmatrixBlockMap(LinearMap.(_get_params_smatrix(sheet, modes, pol))...)
end

function _1DsheetsmatrixBlockMap(halfσˣˣ::LinearMap, halfσʸʸ::LinearMap, pol_kz::LinearMap)::BlockMap
    N = size(pol_kz, 1)
    _fft = LinearMap(x -> fft(x), N)
    _ifft = LinearMap(x -> ifft(x), N)
    A = _get_BlockMap_scattered_side(_fft, _ifft, halfσˣˣ, halfσʸʸ, pol_kz)
    # A = _get_LinearMap_scattered_side(halfσˣˣ, halfσʸʸ, pol_kz)
    B = _get_BlockMap_incident_side(_fft, _ifft, halfσˣˣ, halfσʸʸ, pol_kz)
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

function _get_BlockMap_scattered_side(_fft, _ifft, halfσˣˣ, halfσʸʸ, pol_kz)
    [
        I + _ifft * halfσˣˣ * _fft * pol_kz     -I - _ifft * halfσˣˣ * _fft * pol_kz;
        pol_kz + _ifft * halfσʸʸ * _fft         pol_kz + _ifft * halfσʸʸ * _fft
    ]
end

function _get_BlockMap_incident_side(_fft, _ifft, halfσˣˣ, halfσʸʸ, pol_kz)
    [
        -I + _ifft * halfσˣˣ * _fft * pol_kz    I + _ifft * halfσˣˣ * _fft * pol_kz;
        pol_kz - _ifft * halfσʸʸ * _fft         pol_kz - _ifft * halfσʸʸ * _fft
    ]
end