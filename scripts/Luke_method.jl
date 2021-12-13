import LinearAlgebra: I, Diagonal
import FFTW: fftfreq, fft, ifft

"""
    compute_RCWA_method(η::Function, μ::Function, L::Float64, N::Int, k::Float64, i::Int)

Compute the fields using Luke's RCWA method.
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
x:: the positions of the grid points
kz₁:: the values of the wavevector in the z component for y<0
kz₂:: the values of the wavevector in the z component for y>0
Ĩ₁:: an array of the incident fields in the Fourier basis for y<0
Õ₁:: an array of the outgoing fields in the Fourier basis for y<0
Ĩ₂:: an array of the incident fields in the Fourier basis for y>0
Õ₂:: an array of the outgoing fields in the Fourier basis for y>0
"""
function compute_RCWA_method(η::Function, μ::Function, L::Float64, N::Int, k::Float64, i::Int)
    # discretization parameters
    dx = L/N
    x⃗ = range(0, step=dx, length=N)
    k⃗ = 2*pi*fftfreq(N, 1/dx)
    β⃗ = @. sqrt(Complex(k^2 - k⃗^2))
    # construct Fourier-basis operators
    β̂ = Diagonal(β⃗)
    η̃ = ifft(fft(Matrix(Diagonal(η.(x⃗))), 2), 1)
    μ̃ = ifft(fft(Matrix(Diagonal(μ.(x⃗))), 2), 1)
    # build GSTC matrices and scattering matrix
    A = [
        -β̂-k*η̃    -β̂-k*η̃;
        -β̂-k*μ̃     β̂+k*μ̃
    ]
    B = [
        -β̂+k*η̃    -β̂+k*η̃;
        -β̂+k*μ̃     β̂-k*μ̃
    ]
    S = A\B
    # solve 
    Ĩ₁ = [n == i ? 1 : 0 for n in 1:N]
    Ĩ₂ = zeros(N)
    inc = [Ĩ₁; Ĩ₂]
    scat = S*inc
    Õ₁, Õ₂ = scat[1:N], scat[N+1:2N]
    (x=x⃗, kz₁=β⃗, kz₂=β⃗, Ĩ₁=Ĩ₁, Õ₁=Õ₁, Ĩ₂=Ĩ₂, Õ₂=Õ₂)
end