import DeltaRCWA

struct FunctionSheet{A <: Function, B <: Function} <: DeltaRCWA.RCWASheet
    η::A
    μ::B
end

DeltaRCWA.Zₑˣˣ(sheet::FunctionSheet, x) = sheet.μ(x)
DeltaRCWA.Yₘʸʸ(sheet::FunctionSheet, x) = sheet.η(x)

"""
    compute_DeltaRCWA_method(η::Function, μ::Function, L::Float64, N::Int, k::Float64, i::Int)

Compute the fields using DeltaRCWA (equivalent to Luke's RCWA method).
The output fields are on a grid compatible with the BIE method.

Arguments:
η:: A periodic function returning the metasurface magnetic conductance
μ:: A periodic function returning the metasurface electric resistance
L:: the length of the unit cell (some multiple of the period of η, μ)
N:: the number of discretization points along the unit cell
k:: the wavenumber/frequency of the incident planewave (assumed TM polarization)
i:: the index of the incident planewave in the Fourier basis

Returns:
NamedTuple with fields:
I₁:: an array of the incident fields in the position basis for y<0
O₁:: an array of the outgoing fields in the position basis for y<0
I₂:: an array of the incident fields in the position basis for y>0
O₂:: an array of the outgoing fields in the position basis for y>0
"""
function compute_DeltaRCWA_method(η::Function, μ::Function, L::Float64, N::Int, k::Float64, i::Int)
    prob = DeltaRCWA.DeltaRCWAProblem(
        FunctionSheet(η, μ),
        ((N, L), ),
        k,
        DeltaRCWA.TM(),
        [i == j ? 1.0 : 0.0 for j in 1:N],
        zeros(N)
    )
    sol = DeltaRCWA.solve(prob)
    y⃗ = range(-L, L, length=2N)
    y⃗₁ = y⃗[y⃗ .< 0]
    y⃗₂ = y⃗[y⃗ .> 0]
    β⃗ = prob.modes.kz
    O₁ = rotr90(bfft(exp.(-β⃗ * transpose(im * y⃗₁)) .* sol.O₁, 1))
    I₂ = rotr90(bfft(exp.(-β⃗ * transpose(im * y⃗₂)) .* sol.I₂, 1))
    I₁ = rotr90(bfft(exp.( β⃗ * transpose(im * y⃗₁)) .* sol.I₁, 1))
    O₂ = rotr90(bfft(exp.( β⃗ * transpose(im * y⃗₂)) .* sol.O₂, 1))
    (I₁=I₁, O₁=O₁, I₂=I₂, O₂=O₂)
end