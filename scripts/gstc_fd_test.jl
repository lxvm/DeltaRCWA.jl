θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction
θᵗ = π/3 # desired transmitted field angle
d = cos(θᵗ)-cos(θⁱ) # design parameter
k₀ = 10.0 # wavenumber that attains the metasurface design goals
λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals
L₀ = λ₀/abs(d) # period of metasurface
sheet = ComplexExpSheet(θⁱ, k₀, d) # sheet to scatter
L=L₀ # length of domain (some integer multiple of period of metasurface)
k=1.1k₀ # wavenumber of incident field

n = 5
f1(x) = 1.0
f2(x) = 0.0
test_gstc_1d(sheet, k, n, L, f1, f2)