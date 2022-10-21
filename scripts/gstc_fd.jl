#=
Here we use finite differences to test whether the gstcs hold in real space
=#

using DeltaRCWA


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

"""
Test whether the gstcs are satisfied for TM mode scattering in a y-invariant
system
```math
\\begin{align}
-(H_y^{(2,+)} + H_y^{(2,-)} - H_y^{(1,+)} - H_y^{(1,-)})
	&= \\sigma^e_{xx} (E_x^{(2,+)} + E_x^{(2,-)} + E_x^{(1,+)} + E_x^{(1,-)})/2
\\\\
(E_x^{(2,+)} + E_x^{(2,-)} - E_x^{(1,+)} - E_x^{(1,-)})
	&= -\\sigma^m_{yy} (H_y^{(2,+)} + H_y^{(2,-)} + H_y^{(1,+)} + H_y^{(1,-)})/2
\\end{align}
```
where
```math
\\begin{align}
E_x^i &= \\frac{k_z^i}{\\omega \\epsilon^i} H_y^i = \\frac{1}{\\omega \\epsilon^i} \\partial_z H_y^i
\\end{align}
```
will be approximated by finite differences
Input the sheet and number of Fourier modes, and domain size
"""
function test_gstc_1d(sheet, ω, n, L, f1, f2)
    dims = (n, L)
    pw = Planewaves(ω, dims)
    I₁ = DFT * f1.(pw.x⃗[1])
    I₂ = DFT * f2.(pw.x⃗[1])
    sol = solve(DeltaRCWAProblem(sheet, dims, I₁, I₂))
    @assert I₁ ≈ iDFT * sol.I₁
    @assert I₂ ≈ iDFT * sol.I₂
    O₁ = iDFT * sol.O₁
    O₂ = iDFT * sol.O₂
end