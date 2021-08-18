export smatrix

# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers

"""
    smatrix(::RCWAScatterer, modes, ::AbstractPolarization)

Returns the scattering matrix of a sheet acting on the propagating plane wave modes.
The default method returns the identity scattering matrix.
"""
function smatrix(::RCWAScatterer, modes, pol)
    # total number of propagating modes
    # N = sum(modes.is_propagating)
    N = length(modes.kz)
    mortar(
        ((false*I)(N),  I(N)),
        (I(N),          (false*I)(N)),
    )
end

"""
    smatrix(sheet::RCWASheet{1}, modes, pol::UncoupledPolarization)

Returns a 2x2 BlockMatrix for the scattering of modes
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
    # tilde implies Fourier basis
    halfσ̃ˣˣ = 0.5ifft(fft(Matrix(Diagonal(reshape(σˣˣ, :))), 2), 1)
    halfσ̃ʸʸ = 0.5ifft(fft(Matrix(Diagonal(reshape(σʸʸ, :))), 2), 1)
    pol_kz = Diagonal(reshape(modes.kz, :)) / pol_ω
    halfσ̃ˣˣ, halfσ̃ʸʸ, pol_kz
end

"construct the dense scattering matrix"
function _1Dsheetsmatrix(halfσ̃ˣˣ, halfσ̃ʸʸ, pol_kz)
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