import DeltaRCWA

struct FunctionSheet{A <: Function, B <: Function} <: DeltaRCWA.Sheet
    η::A
    μ::B
end

DeltaRCWA.Zₑˣˣ(sheet::FunctionSheet, x) = 0.5sheet.μ(x[1])
DeltaRCWA.Yₘʸʸ(sheet::FunctionSheet, x) = 2sheet.η(x[1])

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
x:: the positions of the grid points
kz₁:: the values of the wavevector in the z component for y<0
kz₂:: the values of the wavevector in the z component for y>0
Ĩ₁:: an array of the incident fields in the Fourier basis for y<0
Õ₁:: an array of the outgoing fields in the Fourier basis for y<0
Ĩ₂:: an array of the incident fields in the Fourier basis for y>0
Õ₂:: an array of the outgoing fields in the Fourier basis for y>0
"""
function compute_DeltaRCWA_method(η::Function, μ::Function, L::Float64, N::Int, k::Float64, i::Int)
    prob = DeltaRCWA.DeltaRCWAProblem(
        FunctionSheet(η, μ),
        ((N, L), ),
        k,
        [i == j ? 1.0 : 0.0 for j in 1:N],
        zeros(N),
    )
    sol = DeltaRCWA.solve(prob)
    x = sol.pw.x⃗[1]
    kz₁ = DeltaRCWA._get_kz(sol.pw, sol.stack.media[1])
    kz₂ = DeltaRCWA._get_kz(sol.pw, sol.stack.media[end])
    (x=x, kz₁=kz₁, kz₂=kz₂, Ĩ₁=sol.I₁, Õ₁=sol.O₁, Ĩ₂=sol.I₂, Õ₂=sol.O₂)
end