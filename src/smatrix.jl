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
    ω, σˣˣ, σʸʸ = _getσωbypol(sheet, modes, pol)
    # tilde implies Fourier basis
    σ̃ˣˣ = ifft(fft(Matrix(Diagonal(reshape(σˣˣ, :))), 2), 1)
    σ̃ʸʸ = ifft(fft(Matrix(Diagonal(reshape(σʸʸ, :))), 2), 1)
    _1Dsheetsmatrix(Diagonal(reshape(modes.kz, :)), ω, σ̃ˣˣ, σ̃ʸʸ)
end

"convenience function to solve problems for each polarization"
_getσωbypol(sheet, modes, ::TE) = (modes.ω * modes.M.μ, σₘˣˣ(sheet, modes.x⃗), σₑʸʸ(sheet, modes.x⃗))
_getσωbypol(sheet, modes, ::TM) = (modes.ω * modes.M.ϵ, σₑˣˣ(sheet, modes.x⃗), σₘʸʸ(sheet, modes.x⃗))

"construct the dense scattering matrix"
function _1Dsheetsmatrix(kz, ω, σ̃ˣˣ, σ̃ʸʸ)
    mortar(
        (-(-I + (σ̃ˣˣ/2) * (kz/ω)),  -I + (σ̃ˣˣ/2) * (kz/ω)),
        ((kz/ω) - (σ̃ʸʸ/2),          (kz/ω) - (σ̃ʸʸ/2)),
    ) \ mortar(
        (-(I + (σ̃ˣˣ/2) * (kz/ω)),   I + (σ̃ˣˣ/2) * (kz/ω)),
        ((kz/ω) + (σ̃ʸʸ/2),          (kz/ω) + (σ̃ʸʸ/2)),
    )
end