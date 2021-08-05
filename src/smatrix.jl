export smatrix

# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers

"""
    smatrix(::T where T <: RCWAScatterer, modes::PlanewaveModes)

Returns the scattering matrix of a sheet acting on the propagating plane wave modes.
The default method returns the identity scattering matrix.
"""
function smatrix(::T where T <: RCWAScatterer, modes::PlanewaveModes, ::AbstractPolarization)::BlockMatrix
    # total number of propagating modes
    # N = sum(modes.is_propagating)
    N = length(modes.kz)
    mortar(
        ((false*I)(N),  I(N)),
        (I(N),          (false*I)(N)),
    )
end

"""
    smatrix(sheet::RCWASheet{1}, modes::PlanewaveModes, pol::UncoupledPolarization)

Returns a 2x2 BlockMatrix for the scattering of modes
"""
function smatrix(sheet::T where T <: RCWASheet{1}, modes::PlanewaveModes, ::TE)
    σˣˣ = ifft(fft(Diagonal(reshape(σₘˣˣ(sheet, modes.x⃗), :)), 2), 1)
    σʸʸ = ifft(fft(Diagonal(reshape(σₑʸʸ(sheet, modes.x⃗), :)), 2), 1)
    _1Dsheetsmatrix(modes.ω * modes.M.μ, reshape(modes.kz, :), σˣˣ, σʸʸ)
end

function smatrix(sheet::T where T <: RCWASheet{1}, modes::PlanewaveModes, ::TM)
    σˣˣ = ifft(fft(Diagonal(reshape(σₑˣˣ(sheet, modes.x⃗), :)), 2), 1)
    σʸʸ = ifft(fft(Diagonal(reshape(σₘʸʸ(sheet, modes.x⃗), :)), 2), 1)
    _1Dsheetsmatrix(modes.ω * modes.M.ϵ, reshape(modes.kz, :), σˣˣ, σʸʸ)
end

function _1Dsheetsmatrix(ω, kz, σˣˣ, σʸʸ)
    mortar(
        (-I + σˣˣ * Diagonal(kz) / (2ω),    I - σˣˣ * Diagonal(kz) / (2ω)),
        (-Diagonal(kz) / ω + σʸʸ/2,         -Diagonal(kz) / ω + σʸʸ/2),
    ) \ mortar(
        (I + σˣˣ * Diagonal(kz) / (2ω),     -(I + σˣˣ * Diagonal(kz) / (2ω))),
        (-(Diagonal(kz) / ω + σʸʸ/2),       -(Diagonal(kz) / ω + σʸʸ/2)),
    )
end