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

L = 1.0
pw = PlaneWaves(20.0, ((30, L),))

# create a lossless surface impedance (purely imaginary)
admitt = rand(ComplexF64, 10)
admitt .= admitt .- conj.(admitt)
impede = rand(ComplexF64, 10)
impede .= impede .- conj.(impede)
sheet = PixelSheet(admitt, impede, (L,))

prob = DeltaRCWA.ScatteringProblem(sheet, pw)
sol = solve(prob)

# retrieve the unitary scattering matrix of between propagating modes in ports
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
    abs(t12)
end

using Zygote
Zygote.gradient(paramloss, admitt, impede)