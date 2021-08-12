export smatrix

# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers

"""
    smatrix(::T where T <: RCWAScatterer, modes::PlanewaveModes)

Returns the scattering matrix of a sheet acting on the propagating plane wave modes.
The default method returns the identity scattering matrix.
"""
function smatrix(::RCWAScatterer, modes::PlanewaveModes, ::AbstractPolarization)::BlockMatrix
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
function smatrix(sheet::RCWASheet{1}, modes::PlanewaveModes, ::TE)
    # tilde implies Fourier basis
    σ̃ˣˣ = ifft(fft(Matrix(Diagonal(reshape(σₘˣˣ(sheet, modes.x⃗), :))), 2), 1)
    σ̃ʸʸ = ifft(fft(Matrix(Diagonal(reshape(σₑʸʸ(sheet, modes.x⃗), :))), 2), 1)
    _1Dsheetsmatrix(modes.ω * modes.M.μ, Diagonal(reshape(modes.kz, :)), σ̃ˣˣ, σ̃ʸʸ)
end

function smatrix(sheet::RCWASheet{1}, modes::PlanewaveModes, ::TM)
    σ̃ˣˣ = ifft(fft(Matrix(Diagonal(reshape(σₑˣˣ(sheet, modes.x⃗), :))), 2), 1)
    σ̃ʸʸ = ifft(fft(Matrix(Diagonal(reshape(σₘʸʸ(sheet, modes.x⃗), :))), 2), 1)
    _1Dsheetsmatrix(modes.ω * modes.M.ϵ, Diagonal(reshape(modes.kz, :)), σ̃ˣˣ, σ̃ʸʸ)
end

function _1Dsheetsmatrix(
    ω::Float64,
    kz::Diagonal{ComplexF64, Array{ComplexF64, 1}},
    σ̃ˣˣ::Array{ComplexF64, 2},
    σ̃ʸʸ::Array{ComplexF64, 2},
)
    mortar(
        (-(-I + (σ̃ˣˣ/2) * (kz/ω)),  -I + (σ̃ˣˣ/2) * (kz/ω)),
        ((kz/ω) - (σ̃ʸʸ/2),          (kz/ω) - (σ̃ʸʸ/2)),
    ) \ mortar(
        (-(I + (σ̃ˣˣ/2) * (kz/ω)),   I + (σ̃ˣˣ/2) * (kz/ω)),
        ((kz/ω) + (σ̃ʸʸ/2),          (kz/ω) + (σ̃ʸʸ/2)),
    )
end