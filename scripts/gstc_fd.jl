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

function DFT(n)
    ω = cispi(-2/n)
    roots = Vector{typeof(ω)}(undef, n)
    roots[1] = one(ω)
    for i in 2:n
        roots[i] = ω*roots[i-1]
    end
    [roots[mod(k*j, n)+1] for k in 0:(n-1), j in 0:(n-1)]
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
The functions f1(x) and f2(x) should compute the incident complex Hy fields at both
sides of the surface in the position basis. They should be periodic with period L
"""
function test_gstc_1d(sheet, ω, n, L, f1, f2)
    dims = (n, L)
    dft = DFT(n)
    pw = Planewaves(ω, dims)
    Hy₁_inc = dft * f1.(pw.x⃗[1])
    Hy₂_inc = dft * f2.(pw.x⃗[1])
    sol = solve(DeltaRCWAProblem(sheet, dims, Hy₁_inc, Hy₂_inc))
    @assert I₁ ≈ idft * sol.I₁
    @assert I₂ ≈ idft * sol.I₂
    Hy₁_out = (sol.O₁' * dft)' ./ n # idft * sol.O₁
    Hy₂_out = (sol.O₂' * dft)' ./ n # idft * sol.O₂
end

"""
Use finite differences to estimate
```math
\\begin{align}
E_x^i &= \\frac{k_z^i}{\\omega \\epsilon^i} H_y^i = \\frac{1}{\\sqrt{-1} \\omega \\epsilon^i} \\partial_z H_y^i
\\end{align}
```
Specifically, uses a 1st order, one-sided 
"""
function fd(Hy, L)
end