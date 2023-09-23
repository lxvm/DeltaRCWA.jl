using DeltaRCWA
using LinearAlgebra
using Statistics

struct PixelSheet{N,T} <: Sheet
    admitt::Array{Complex{T},N}
    impede::Array{Complex{T},N}
    L::NTuple{N,T}
end

function DeltaRCWA.Zₑˣˣ(sheet::PixelSheet{N}, x⃗::NTuple{N}) where {N}
    idx = map((x,l,s) -> mod(floor(Int, s*x/l), s)+1, x⃗, sheet.L, size(sheet.impede))
    sheet.impede[idx...]
end

function DeltaRCWA.Yₘʸʸ(sheet::PixelSheet{N}, x⃗::NTuple{N}) where {N}
    idx = map((x,l,s) -> mod(floor(Int, s*x/l), s)+1, x⃗, sheet.L, size(sheet.admitt))
    sheet.admitt[idx...]
end

L = 1.0     # unit cell period
nmode = 30  # number of Fourier modes in frequency discretization
ω = 2.0     # frequency
pw = PlaneWaves(ω, ((nmode, L),))

# create a lossless surface impedance (purely imaginary)
npix = 10   # pixels to use in impedance sheet
admitt = rand(ComplexF64, npix)
admitt .= admitt .- conj.(admitt)
impede = rand(ComplexF64, npix)
impede .= impede .- conj.(impede)
sheet = PixelSheet(admitt, impede, (L,))

prob = DeltaRCWA.ScatteringProblem(sheet, pw)
sol = solve(prob)

# retrieve the unitary scattering matrix of propagating modes across ports
U, = DeltaRCWA.propagating_unitary_smatrix(sol)
# norm(U'U - I)/ norm(U'U) ≈ 1 # for lossless surface impedances

# retrieve the normalized transmission coefficients between normal-incidence
# modes of each port
t12, t21 = DeltaRCWA.transmission_coefficient(sol)


function paramloss(admitt, impede)
    L = 1.0
    sheet = PixelSheet(admitt, impede, (L,))
    pw = PlaneWaves(20.0, ((30, L),))
    prob = DeltaRCWA.ScatteringProblem(sheet, pw)
    sol = solve(prob)
    t12, t21 = DeltaRCWA.transmission_coefficient(sol)
    abs(t12)    # scalar function of transmission
end

using Zygote
Zygote.gradient(paramloss, admitt, impede)