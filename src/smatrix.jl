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
    σˣˣ = ifft(fft(Matrix(Diagonal(reshape(σₘˣˣ(sheet, modes.x⃗), :))), 2), 1)
    σʸʸ = ifft(fft(Matrix(Diagonal(reshape(σₑʸʸ(sheet, modes.x⃗), :))), 2), 1)
    _1Dsheetsmatrix(modes.ω * modes.M.μ, Diagonal(reshape(modes.kz, :)), σˣˣ, σʸʸ)
end

function smatrix(sheet::T where T <: RCWASheet{1}, modes::PlanewaveModes, ::TM)
    σˣˣ = ifft(fft(Matrix(Diagonal(reshape(σₑˣˣ(sheet, modes.x⃗), :))), 2), 1)
    σʸʸ = ifft(fft(Matrix(Diagonal(reshape(σₘʸʸ(sheet, modes.x⃗), :))), 2), 1)
    _1Dsheetsmatrix(modes.ω * modes.M.ϵ, Diagonal(reshape(modes.kz, :)), σˣˣ, σʸʸ)
end

function _1Dsheetsmatrix(
    ω::Float64,
    kz::Diagonal{ComplexF64, Array{ComplexF64, 1}},
    σˣˣ::Array{ComplexF64, 2},
    σʸʸ::Array{ComplexF64, 2},
)
    mortar(
        (-(-I + (σˣˣ/2) * (kz/ω)),  -I + (σˣˣ/2) * (kz/ω)),
        ((kz/ω) - (σʸʸ/2),          (kz/ω) - (σʸʸ/2)),
    ) \ mortar(
        (-(I + (σˣˣ/2) * (kz/ω)),   I + (σˣˣ/2) * (kz/ω)),
        ((kz/ω) + (σʸʸ/2),          (kz/ω) + (σʸʸ/2)),
    )
end