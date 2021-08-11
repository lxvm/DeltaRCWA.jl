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
using LinearAlgebra: Diagonal, I

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

# ╔═╡ db5067ad-60a3-4f75-a25e-441ccf61ea6f
md"
## Using DeltaRCWA
In this section we define and solve a scattering problem that identical to Luke's
(except for the difference in convention regarding conductivities vs impedances).
I will also give this material these properties for the TM polarization,
since this would cause infinite electric conductivity but finite magnetic.
These are the stages to using the solver, which needs this data
```julia
struct DeltaRCWAProblem{N}
    structure::RCWAScatterer{N, M} where M
    modes::PlanewaveModes{N}
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
This data can be supplied with the following information:
#### Wave frequency
Choose an `ω::Float64` for the temporal frequency of the wave
#### Unit cell dimensions and discretization
For each periodic dimension, give a number of grid points per unit cell and a size of 
the unit cell: `dims::NTuple{N, Tuple{Int64,Float64}}`
#### Wave medium
Planewaves are the eigenmodes of uniform media. Define the medium properties with the
constructor `UniformMedium{T <: Number}(ϵ::T, μ::T)`, or for vacuum just `Vacuum()`.
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
#### Methods that a `RCWASheet` implements
Any `RCWASheet` may implement the following 8 methods
`σₑˣˣ, σₑˣʸ, σₑʸˣ, σₑʸʸ, σₘˣˣ, σₘˣʸ, σₘʸˣ, σₘʸʸ`
one for each component of the electric and magnetic conductivity matrices.
These methods return the value of that conductivity component at all points in the
unit cell, whose coordinates are given by a product iterator over the argument `x⃗`.
By default/fallback, these methods all call for a trivial scattering sheet
```julia
function nonconducting(::RCWASheet{N}, x⃗::NTuple{N, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}) where N
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
`plot(::DeltaRCWASolution{1})` should just work! See the examples below for usage.
### Project TODO
- Implement a solver for `RCWASheet{2}`/3D photonic crystals 
- Implement functions to:
  - Compute the reflected and transmitted power
  - Compute the complex transmission coefficient
- Fast Redheffer star product
- Benchmarking and testing
- Differentiability of solver
"

# ╔═╡ d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
modeN = 0

# ╔═╡ aa87c7bd-4d9e-49da-8cf0-4a0fb6081588
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
N₀(x) = @. -sin(θ)*(1-exp(1im*k*d*x))

nvec = 0:99
dx = L/length(nvec)
xvec = [n*dx-L/2 for n in nvec]
kₓ = 2*pi*fftfreq(length(nvec), 1/dx)
β = @. sqrt(Complex(k^2 - kₓ^2))

u_p = ComplexF64[n == modeN ? 1 : 0 for n in nvec]
# u_p = fft(uInc(xvec, 0))/length(nvec)
u_n = zeros(ComplexF64, length(nvec))
M = Matrix(Diagonal(M₀(xvec)))
N = Matrix(Diagonal(N₀(xvec)))

A = [-Diagonal(β)-k*ifft(fft(M, 2), 1)    -Diagonal(β)-k*ifft(fft(M, 2), 1);
    -Diagonal(β)-k*ifft(fft(N, 2), 1)     Diagonal(β)+k*ifft(fft(N, 2), 1)]
B = [-Diagonal(β)+k*ifft(fft(M, 2), 1)    -Diagonal(β)+k*ifft(fft(M, 2), 1);
    -Diagonal(β)+k*ifft(fft(N, 2), 1)     Diagonal(β)-k*ifft(fft(N, 2), 1)]
S = A\B

u_in = [u_p; u_n]
u_out = S*u_in
u_out_p = u_out[1:length(nvec)]
u_out_n = u_out[length(nvec)+1:2*length(nvec)]

E1(x,z) = sum([u_out_p[n+1]*exp(-1im*kₓ[n+1]*x)*exp(1im*β[n+1]*z) for n in nvec])+sum([u_p[n+1]*exp(-1im*kₓ[n+1]*x)*exp(-1im*β[n+1]*z) for n in nvec])
E2(x,z) = sum([u_out_n[n+1]*exp(-1im*kₓ[n+1]*x)*exp(-1im*β[n+1]*z) for n in nvec])+sum([u_n[n+1]*exp(-1im*kₓ[n+1]*x)*exp(1im*β[n+1]*z) for n in nvec])

xview = -L:0.02:L 
zview = 0:0.01:L
zview2 = -L:0.01:0
Emat = zeros(Complex{Float64}, 0)
Emat2 = zeros(Complex{Float64}, 0)
for (z1, z2) in zip(zview, zview2)
    append!(Emat, E1.(xview, z1))
    append!(Emat2, E2.(xview, z2))
end
Emat = reshape(Emat,length(xview),:)
Emat2 = reshape(Emat2,length(xview),:)
	
plt1 = heatmap(xview, zview2, transpose(real.(Emat2)))
plt1 = heatmap!(xview, zview, transpose(real.(Emat)))
plt2 = heatmap(xview, zview2, transpose(imag.(Emat2)))
plt2 = heatmap!(xview, zview, transpose(imag.(Emat)))
plt3 = heatmap(xview, zview2, transpose(abs2.(Emat2)))
plt3 = heatmap!(xview, zview, transpose(abs2.(Emat)))
end;

# ╔═╡ c34b5661-e650-4fc1-9dcf-285c39ddb983
plt1 # real
# plt2 # imag
# plt3 # abs2

# ╔═╡ 65231a8c-f54d-4ffc-bf59-9dc4cc33f61a
begin
### DeltaRCWA solver
struct ComplexExpSheet{T <: Real} <: RCWASheet{1}
    θ::T
    θᵗ::T
    k::T
    d::T
    L::T
end

function DeltaRCWA.σₑˣˣ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
    2 ./ ComplexF64[e ≈ 0 ? 1e-14 : e for e in N₀(x⃗...)]
end

function DeltaRCWA.σₘʸʸ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
    2M₀(x⃗...)
end

ω = k
sheet = ComplexExpSheet(θ, θᵗ, k, d, L)
modes = PlanewaveModes(ω, ((length(nvec), sheet.L), ), Vacuum())
pol = TM()
prob = DeltaRCWAProblem(sheet, modes, pol, u_p, u_n)
sol = solve(prob)
end

# ╔═╡ 65cb38cd-89bc-4254-9ebd-c4d583d8cfb7
plot(sol; part=real, method=:fft, combine=false)

# ╔═╡ 786a6947-7082-4902-a625-8be4bd3e30e7
begin
	### Compare the resitivities/conductivities
	plot(xvec,  real.(M₀(modes.x⃗...)), label="Re(M₀)")
	plot!(xvec, imag.(M₀(modes.x⃗...)), label="Im(M₀)")
	plot!(xvec, real.(σₘʸʸ(sheet, modes.x⃗)), label="Re(σₘʸʸ)")
	plot!(xvec, imag.(σₘʸʸ(sheet, modes.x⃗)), label="Im(σₘʸʸ)")
end

# ╔═╡ 3f69f8c5-4d79-47f4-973e-3286cc9b4f6a
begin
	plot(xvec, real.(N₀(modes.x⃗...)), label="Re(N₀)")
	plot!(xvec, imag.(N₀(modes.x⃗...)), label="Im(N₀)")
	plot!(xvec, real.(σₑˣˣ(sheet, modes.x⃗)), label="Re(σₑˣˣ)")
	plot!(xvec, imag.(σₑˣˣ(sheet, modes.x⃗)), label="Im(σₑˣˣ)", ylims=(-5, 5))
end

# ╔═╡ a260a293-bcb4-483c-ba20-2af3e3cbfa58
begin
	# plot(real.(sol.I₁))
	# plot(real.(Luke_sol.I₁))
	# plot(real.(sol.O₂))
	# plot(real.(ifft(sol.O₂)))
	# plot(real.(Luke_sol.O₂))
	# plot(real.(ifft(Luke_sol.O₂)))
	# plot(real.(Luke_sol.O₁))
	# plot(real.(ifft(Luke_sol.O₁)))
end

# ╔═╡ 3eb5ee6e-8581-41d5-9b72-33d82bad4c5b
# Check that the outputs and inputs match
u_out_p ≈ sol.O₂

# ╔═╡ fc55344b-f0ff-44da-812c-43255aa4e7fd
u_out_n ≈ sol.O₁

# ╔═╡ 51f5b8a1-f5ec-4b42-a2b3-c3bb9f33e063
u_p ≈ sol.I₁

# ╔═╡ 36afdf48-e275-451b-a095-f2e5e57e5167
u_n ≈ sol.I₂

# ╔═╡ 903b1a7b-8f92-451f-92d7-7e8a72a8c7aa
md"
Here we need to map Luke's convention to mine.
- Luke's choice of side/sign is opposite of mine
  - I chose the +z direction to go from port 1 to port 2 but Luke is opposite
  - the _p and _n subscripts are consistent with port 1 and port 2 according to Luke's derivations
"

# ╔═╡ 2d1f452b-4d01-4fca-ae65-a864c4afa842
Luke_sol = DeltaRCWASolution(modes, pol, u_p, u_n, u_out_p, u_out_n)

# ╔═╡ a10c8dee-6db7-475e-89f2-bc3c064522f9
plot(Luke_sol; part=real, method=:fft, combine=false)

# ╔═╡ a9b12d46-c30d-4fca-9ce5-d3dd7498cf5d
# check that the S-matrix matches (with or without blockarray)
S ≈ smatrix(sheet, modes, pol)

# ╔═╡ ca68a6e9-6697-43c3-a51c-966c17b101cc
# plot the spectrum of the matrix
plot(S)

# ╔═╡ d676e456-9984-4d24-8c8f-fece0edb1ee8
plot(Matrix(smatrix(sheet, modes, pol)))

# ╔═╡ Cell order:
# ╠═93f34ece-0216-4722-9bdc-70ee684d9bd3
# ╠═628d5d5d-2753-47e3-a02b-21b4a89a159e
# ╠═9c58eb1a-313f-4f89-9958-7f33c7a072d3
# ╠═b4a594f1-a38c-4ede-9860-d4a5afae15c5
# ╠═55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
# ╠═2a30d4a7-ef57-478d-93ac-b5756b6f3909
# ╠═adbd86c3-f970-4681-bcde-ddda1050eefd
# ╠═67fb8117-7d5b-4536-9e36-7dda36997dff
# ╠═d915d60d-b159-400a-811c-af9b8828ec91
# ╠═aa87c7bd-4d9e-49da-8cf0-4a0fb6081588
# ╠═c34b5661-e650-4fc1-9dcf-285c39ddb983
# ╟─db5067ad-60a3-4f75-a25e-441ccf61ea6f
# ╠═65231a8c-f54d-4ffc-bf59-9dc4cc33f61a
# ╠═65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
# ╠═786a6947-7082-4902-a625-8be4bd3e30e7
# ╠═3f69f8c5-4d79-47f4-973e-3286cc9b4f6a
# ╠═65cb38cd-89bc-4254-9ebd-c4d583d8cfb7
# ╠═d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
# ╠═a260a293-bcb4-483c-ba20-2af3e3cbfa58
# ╠═3eb5ee6e-8581-41d5-9b72-33d82bad4c5b
# ╠═fc55344b-f0ff-44da-812c-43255aa4e7fd
# ╠═51f5b8a1-f5ec-4b42-a2b3-c3bb9f33e063
# ╠═36afdf48-e275-451b-a095-f2e5e57e5167
# ╠═903b1a7b-8f92-451f-92d7-7e8a72a8c7aa
# ╠═2d1f452b-4d01-4fca-ae65-a864c4afa842
# ╠═a10c8dee-6db7-475e-89f2-bc3c064522f9
# ╠═a9b12d46-c30d-4fca-9ce5-d3dd7498cf5d
# ╠═ca68a6e9-6697-43c3-a51c-966c17b101cc
# ╠═d676e456-9984-4d24-8c8f-fece0edb1ee8
