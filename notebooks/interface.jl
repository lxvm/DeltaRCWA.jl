### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 628d5d5d-2753-47e3-a02b-21b4a89a159e
import Pkg

# ╔═╡ 9c58eb1a-313f-4f89-9958-7f33c7a072d3
Pkg.activate(".")

# ╔═╡ b4a594f1-a38c-4ede-9860-d4a5afae15c5
using Revise

# ╔═╡ 55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
using DeltaRCWA

# ╔═╡ 2a30d4a7-ef57-478d-93ac-b5756b6f3909
using LinearAlgebra: Diagonal, I, norm

# ╔═╡ 805f9d8b-ae21-4b27-9d6d-6a818b1fabc9
using BlockArrays: BlockMatrix

# ╔═╡ adbd86c3-f970-4681-bcde-ddda1050eefd
using FFTW: fftfreq, fft, ifft

# ╔═╡ 67fb8117-7d5b-4536-9e36-7dda36997dff
using BlockArrays: Block, mortar

# ╔═╡ 65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
using Plots

# ╔═╡ 93f34ece-0216-4722-9bdc-70ee684d9bd3
md"
# Interface
This notebook shows how to use the DeltaRCWA interface for 2D photonic crystals
and compares the output of the package to a direct calculation as done by Luke.
"

# ╔═╡ d915d60d-b159-400a-811c-af9b8828ec91
md"
## Direct calculation
This is a copy of Luke's code testing Carlos' surface impedance example, including
the fields after scattering.
"

# ╔═╡ e9534450-9a3d-4efa-a3e6-4c9aba4e1646
begin
#### Luke's solver
k = 10.0
λ = 2*pi/k
θ = -π/2.0 # incidence angle (measured with respect to the x axis)
α = k*cos(θ)
β = k*sin(θ)
uInc(x,y)= @. exp(1im*α*x+1im*β*y)  # incident planewave

θᵗ = -π/8    # transmitted field angle
d  = cos(θᵗ)-cos(θ)
L = 2*(2*π)/(k*abs(d))  # Unit cell width
M₀(x) = @. -sin(θ)*(1+exp(1im*k*d*x))
N₀(x) = @. 1e-8-sin(θ)*(1-exp(1im*k*d*x))

nvec = 0:99
dx = L/length(nvec)
xvec = [n*dx-L/2 for n in nvec]
kₓ = 2*pi*fftfreq(length(nvec), 1/dx)
β = @. sqrt(Complex(k^2 - kₓ^2))

M = Matrix(Diagonal(M₀(xvec)))
N = Matrix(Diagonal(N₀(xvec)))

A = [-Diagonal(β)-k*ifft(fft(M, 2), 1)    -Diagonal(β)-k*ifft(fft(M, 2), 1);
    -Diagonal(β)-k*ifft(fft(N, 2), 1)     Diagonal(β)+k*ifft(fft(N, 2), 1)]
B = [-Diagonal(β)+k*ifft(fft(M, 2), 1)    -Diagonal(β)+k*ifft(fft(M, 2), 1);
    -Diagonal(β)+k*ifft(fft(N, 2), 1)     Diagonal(β)-k*ifft(fft(N, 2), 1)]

S = A\B
end;

# ╔═╡ ad9136c8-ec61-4736-925a-b4a2165080c6
# Is the scattering matrix unitary?
norm(S'S-I)/norm(S'S)

# ╔═╡ d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
### choose which mode to scatter
modeN = 0

# ╔═╡ db5067ad-60a3-4f75-a25e-441ccf61ea6f
md"
## Using DeltaRCWA
In this section we define and solve a scattering problem that identical to Luke's
(except for the difference in convention regarding conductivities vs impedances).
I will also give this material these properties for the TM polarization,
since this would cause infinite electric conductivity but finite magnetic.
These are the stages to using the solver, which needs this data
```julia
struct DeltaRCWAProblem{T₁, T₂, N, L}
    structure::SheetStack{T₁, N, L}
    modes::PlanewaveModes{T₂, N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
end
```
as explained and demonstrated below:
### Defining the scattering modes
A struct called `PlanewaveModes{N}` is used to specify the discretization and
periodicities of space in the unit cell where `DeltaRCWA` solves for the fields.
`N` refers to the number of dimensions along which the unit cell is periodic.
Along these periodic dimensions, the solutions can be expanded in the Fourier basis
due to Bloch's theorem (much like the [Kronig-Penny model](https://en.wikipedia.org/wiki/Particle_in_a_one-dimensional_lattice)).
The full structure is:
```julia
struct PlanewaveModes{T, N}
    ω::Float64
    M::UniformMedium{T}
    dims::NTuple{N, Tuple{Int64, Float64}}
    x⃗::NTuple{N, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}
    k⃗::NTuple{N, Frequencies{Float64}}
    kz::Array{ComplexF64, N}
end
```
This data can be supplied with the following information (using convenience methods):
#### Wave frequency
Choose an `ω::Float64` for the temporal frequency of the wave
#### Unit cell dimensions and discretization
For each periodic dimension, give a number of grid points per unit cell and a size of 
the unit cell: `dims::NTuple{N, Tuple{Int64,Float64}}`
#### Wave medium
Planewaves are the eigenmodes of uniform media. Define the medium properties with
either constructor
```julia
struct UniformMedium{T <: Number}
    ϵ::T
    μ::T
end

Vacuum() = UniformMedium(true, true)
```
#### Polarization
This isn't stored in `PlanewaveModes`, but it is necessary information for 2D
photonic crystals with diagonal conductivity matrices, when the TE and TM
polarizations decouple. Tell the solver whether to solve `TE()` or `TM()` when
`N=1`. `Coupled()` exists for `N=2` (to be implemented).
#### Incident mode amplitudes
Defining a distribution over the incident planewave modes allows for a solution of the
scattering problem. Two distributions must be specified: one at each port of the device.
### Defining the scattering structure
The abstract type `RCWAScatterer{N,M}` is the main respresentation of a unit cell 
in a periodic scattering structure in `DeltaRCWA`.
The parameters `N,M` are integers to specify that a structure has `N` periodic
dimensions (thus is at least a `N+1`-dimensional photonic crystal, after adding the
scattering axis) and that occupies `M` dimensions.
The useful subtypes of this are `RCWASheet{N} <: RCWAScatterer{N,2}` and
`RCWASlab{N} <: RCWAScatterer{N,3}`, which represent structures with surface
impedances and volume impedances, respectively.
`DeltaRCWA` concerns itself only with the `RCWASheet` since any nontrivial `RCWASlab`
would require a more computationally expensive solver for the eigenmodes of Maxwell's
equations in that structure.
Other software, such as [S⁴](http://www.stanford.edu/group/fan/S4/), can do this more
computationally expensive step, though `DeltaRCWA` could be extended to do it as well
by defining a method for `smatrix(::RCWASlab{N}, ::PlanewaveModes{N}, ...)`.
#### Creating a `RCWASheet`
To create a sheet as part of a 2D photonic crystal, create a struct that is a subtype
of `RCWASheet{1}` and store all the geometric and material parameters you need in your
struct to define the electric and magnetic conductivity matrices for it.
#### Creating a `SheetStack`
In order to put 2 or more sheets in a sequence, separated by a uniform medium
(the one defined in `PlanewaveModes`) create an instance of a `SheetStack`.
Create a Tuple of the sheets you want to scatter off of and create a second Tuple
with the size of the Vacuum gap that separates each of the sheets and pass these to
the `SheetStack` constructor. Note that there is one gap fewer than the number of sheets.
```julia
struct SheetStack{T, N, L} <: RCWAScatterer{T, N}
    sheets::Tuple{RCWASheet{T, N}, Vararg{RCWASheet{T, N}, L}}
    depths::Tuple{Vararg{Float64, L}}
end
```
For convenience in defining a problem with one sheet, the DeltaRCWAProblem constructor
will turn the sheet into a sheetstack for you.
#### Methods that a `RCWASheet` implements
Any `RCWASheet` may implement the following 8 methods
`σₑˣˣ, σₑˣʸ, σₑʸˣ, σₑʸʸ, σₘˣˣ, σₘˣʸ, σₘʸˣ, σₘʸʸ`
one for each component of the electric and magnetic conductivity matrices.
These methods return the value of that conductivity component at all points in the
unit cell, whose coordinates are given by a product iterator over the argument `x⃗`.
By default/fallback, these methods all call for a trivial scattering sheet
```julia
function nonconducting(::RCWASheet, x⃗) where N
    zeros(Bool, length.(x⃗))
end
```
unless you define a method which dispatches on your type, i.e.:
```julia
σₑˣˣ(sheet::MySheet, x⃗) = ...
```
See the example below on how to extend the methods exported by `DeltaRCWA`, 
or read the [documentation](https://docs.julialang.org/en/v1/manual/modules/#using-and-import-with-specific-identifiers,-and-adding-methods)
The output of these methods can only depend on the input point and the geometric
parameters stored in the sheet type.
These methods are used to construct the `smatrix` for each sheet using the formula
defined [here](https://github.com/lxvm/DeltaRCWA.jl/blob/main/src/smatrix.jl).
### Problem and solution
Create the `DeltaRCWAProblem` object, `solve` it, and receive a `DeltaRCWASolution`
object.
### Analysis
The `DeltaRCWASolution` objects have special methods to analyze and visualize them.
#### Plotting
`plot(::DeltaRCWASolution{T, 1} where T)` should just work!
See the examples below for usage.
### Project TODO
- Implement a solver for `RCWASheet{2}`/3D photonic crystals 
- Implement functions to:
  - Compute the reflected and transmitted power
  - Compute the complex transmission coefficient
- Fast Redheffer star product
- Benchmarking and testing
- Differentiability of solver
"

# ╔═╡ 4125d1a4-2a57-431a-b7ea-ab8f44994143
begin
### DeltaRCWA solver
struct ComplexExpSheet{T} <: RCWASheet{T, 1}
    θ::T
    θᵗ::T
    k::T
    d::T
    L::T
end

### Define how to convert between M/N and conductivity matrix conventions
function DeltaRCWA.σₑˣˣ(sheet::ComplexExpSheet, x⃗)
	# 2 ./ N₀(x⃗...)
    2 ./ ComplexF64[e ≈ 0 ? 1e-14 : e for e in N₀(x⃗...)]
end

function DeltaRCWA.σₘʸʸ(sheet::ComplexExpSheet, x⃗)
    2M₀(x⃗...)
end

ω = k
sheet = ComplexExpSheet(θ, θᵗ, k, d, L)
dims = ((length(nvec), sheet.L), )
modes = PlanewaveModes(ω, dims, Vacuum())
pol = TM()
end;

# ╔═╡ 786a6947-7082-4902-a625-8be4bd3e30e7
begin
	### Display the magnetic conductivity (M) / electric resistivity (N) along sheet
	plot(xvec,  real.(M₀(modes.x⃗...)), label="Re(M₀)")
	plot!(xvec, imag.(M₀(modes.x⃗...)), label="Im(M₀)", linestyle=:dash)
	plot!(xvec, real.(N₀(modes.x⃗...)), label="Re(N₀)")
	plot!(xvec, imag.(N₀(modes.x⃗...)), label="Im(N₀)", linestyle=:dash)
end

# ╔═╡ 483e04a7-ac35-44e1-88e7-6e18737d7110
begin
u_p = ComplexF64[n == modeN ? 1 : 0 for n in nvec]
# u_p = fft(uInc(xvec, 0))/length(nvec) # too noisy to plot
u_n = zeros(ComplexF64, length(nvec))
u_in = [u_p; u_n]
u_out = S*u_in
u_out_p = u_out[1:length(nvec)]
u_out_n = u_out[length(nvec)+1:2*length(nvec)]
Luke_sol = DeltaRCWASolution(modes, pol, u_p, u_n, u_out_p, u_out_n)
end;

# ╔═╡ 2d1f452b-4d01-4fca-ae65-a864c4afa842
plot(Luke_sol; part=real, method=:fft, combine=false)

# ╔═╡ 0a3d4a56-d783-4bd9-9392-17ade6242a97
### verify that the scattering matrices match
S ≈ smatrix(sheet, modes, pol)

# ╔═╡ 77400f50-4e25-4fe6-8a0a-16f6cf6cb150
begin
prob = DeltaRCWAProblem(sheet, modes, pol, u_p, u_n)
sol = solve(prob, method=:matrixfree)
plot(sol;)
end

# ╔═╡ 7efe3220-c28b-4162-978a-7cf20673b1c4
begin
### demonstrate SheetStack
struct TrivialSheet{T} <: RCWASheet{T, 1} end
nsheets = 2
gap(x) = 2L
stack = SheetStack(
	Tuple(TrivialSheet{Bool}() for i in 1:nsheets),
	Tuple(gap(i) for i in 1:(nsheets-1)),
)
stackprob = DeltaRCWAProblem(stack, modes, pol, u_p, u_n)
stacksol = solve(stackprob)
end;

# ╔═╡ d8d2123c-92ce-422e-bb68-bb09e689d44c
plot(stacksol; part=imag, combine=true)

# ╔═╡ 418a2244-63dd-4922-93e7-4e34c9cdb583
begin
	### Demonstrate Fresnel scattering
	water = UniformMedium(80.3, 1.0)
	water_modes = PlanewaveModes(ω, dims, water)
	intf = UniformInterface(Vacuum(Float64), water)
	intf_S = smatrix(intf, modes, water_modes)
end;

# ╔═╡ 185390b7-2d55-4010-93eb-c9e33662435b
begin
	intf_out = intf_S * mortar([u_p, u_n])
	int_O₁ = intf_out[1:length(nvec)]
	int_O₂ = intf_out[(1+length(nvec)):end]
	intf_sol = DeltaRCWASolution(modes, pol, u_p, u_n, int_O₁, int_O₂)
end

# ╔═╡ 38b5a713-b224-4d62-baf2-69a1ef93e0bc
plot(intf_sol)

# ╔═╡ Cell order:
# ╠═93f34ece-0216-4722-9bdc-70ee684d9bd3
# ╠═628d5d5d-2753-47e3-a02b-21b4a89a159e
# ╠═9c58eb1a-313f-4f89-9958-7f33c7a072d3
# ╠═b4a594f1-a38c-4ede-9860-d4a5afae15c5
# ╠═55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
# ╠═2a30d4a7-ef57-478d-93ac-b5756b6f3909
# ╠═805f9d8b-ae21-4b27-9d6d-6a818b1fabc9
# ╠═adbd86c3-f970-4681-bcde-ddda1050eefd
# ╠═67fb8117-7d5b-4536-9e36-7dda36997dff
# ╠═65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
# ╠═d915d60d-b159-400a-811c-af9b8828ec91
# ╠═e9534450-9a3d-4efa-a3e6-4c9aba4e1646
# ╠═786a6947-7082-4902-a625-8be4bd3e30e7
# ╠═ad9136c8-ec61-4736-925a-b4a2165080c6
# ╠═d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
# ╠═483e04a7-ac35-44e1-88e7-6e18737d7110
# ╠═2d1f452b-4d01-4fca-ae65-a864c4afa842
# ╠═db5067ad-60a3-4f75-a25e-441ccf61ea6f
# ╠═4125d1a4-2a57-431a-b7ea-ab8f44994143
# ╠═0a3d4a56-d783-4bd9-9392-17ade6242a97
# ╠═77400f50-4e25-4fe6-8a0a-16f6cf6cb150
# ╠═7efe3220-c28b-4162-978a-7cf20673b1c4
# ╠═d8d2123c-92ce-422e-bb68-bb09e689d44c
# ╠═418a2244-63dd-4922-93e7-4e34c9cdb583
# ╠═185390b7-2d55-4010-93eb-c9e33662435b
# ╠═38b5a713-b224-4d62-baf2-69a1ef93e0bc