using DeltaRCWA
using LinearAlgebra
using Statistics

struct ComplexExpSheet{T} <: Sheet
    θ::T # incidence angle
    k::T # wavenumber
    d::T # designer parameter
end

function DeltaRCWA.Zₑˣˣ(sheet::ComplexExpSheet, x⃗::Tuple)
    z = -0.5sin(sheet.θ)*(1-exp(1im*sheet.k*sheet.d*x⃗[1]))
    im*imag(z)
end

function DeltaRCWA.Yₘʸʸ(sheet::ComplexExpSheet, x⃗::Tuple)
    y = -2sin(sheet.θ)*(1+exp(1im*sheet.k*sheet.d*x⃗[1]))
    im*imag(y)
end

θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction
θᵗ = π/3 # desired transmitted field angle
d = cos(θᵗ)-cos(θⁱ) # design parameter
k₀ = 10.0 # wavenumber that attains the metasurface design goals
λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals
L₀ = λ₀/abs(d) # period of metasurface
sheet = ComplexExpSheet(θⁱ, k₀, d) # sheet to scatter
L=L₀ # length of domain (some integer multiple of period of metasurface)
k=1.1k₀ # wavenumber of incident field

for n in range(5, length=10, step=5)
    @show n
    # 2D
    dims = ((n, L), )
    nrep = 2
    # 3D
    # dims = ((n, L), (n,L))
    # nrep = 4
    pw = PlaneWaves(k, dims)
    S = smatrix(HField(), pw, sheet, Vacuum, Vacuum)
    @show size(S)
    kz = DeltaRCWA._get_kz(pw, Vacuum)
    mask = vec(isreal.(kz) .& !iszero(kz))
    @show size(mask) count(mask)
    SS = S[repeat(mask, nrep), repeat(mask, nrep)]
    @show norm(SS'SS - I)/norm(SS'SS)
    ee = eigen(SS)
    @show mean(abs.(ee.values))
    @show var(abs.(ee.values))
    # to Poynting amplitudes
    ϵ₀, μ₀ = DeltaRCWA._get_ϵμ(Vacuum)
    Z₀ = sqrt(μ₀/ϵ₀)
    v₀ = sqrt(1/(μ₀*ϵ₀))
    Y₀ = inv(Z₀)
    @show twiddle = kz[mask]*Z₀/(2*k) #  for HField
    # twiddle = kz[mask]*Y₀/(2*k) #  for EField
    A = Diagonal(repeat(twiddle, nrep))
    SSS = inv(A)*SS*A
    @show norm(SSS'SSS - I)/norm(SSS'SSS)
    @show diag(eigvecs(SS))
end