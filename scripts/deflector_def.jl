const θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction
const θᵗ = π/3 # desired transmitted field angle
const d = cos(θᵗ)-cos(θⁱ) # design parameter
const k₀ = 10.0 # wavenumber that attains the metasurface design goals
const λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals
const L₀ = λ₀/abs(d) # period of metasurface

η(x) = sin(abs(θⁱ))*(1+exp.(im*k₀*d*x))
μ(x) = sin(abs(θⁱ))*(1-exp.(im*k₀*d*x))

const L=L₀ # length of domain (some integer multiple of period of metasurface)
const k=1.1k₀ # wavenumber of incident field

const α = k * cos(θⁱ)
const β = k * sin(θⁱ)
const i = 1 # == 1 for θⁱ = π/2
const exponents = 5:9 # defines the number of grid points, base 2