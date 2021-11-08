### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 628d5d5d-2753-47e3-a02b-21b4a89a159e
import Pkg

# ╔═╡ 9c58eb1a-313f-4f89-9958-7f33c7a072d3
Pkg.activate(".")

# ╔═╡ 55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
using DeltaRCWA

# ╔═╡ 65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
using Plots

# ╔═╡ 0a627414-a5af-4fc3-852f-c98105e4d860
md"
## Using `DeltaRCWA`'s interface
These are the stages to using the solver, which needs this data
```julia
struct DeltaRCWAProblem{N, T₁<:RCWAScatterer{N}, T₂}
    structure::T₁
    modes::PlanewaveModes{T₂, N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
end
```
as explained and demonstrated below:
### Defining the scattering modes
A struct called `PlanewaveModes{T, N}` is used to specify the discretization and
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
"

# ╔═╡ f0ca02cb-c7b0-42d8-b598-98499aa09f12
md"
#### Wave frequency
Choose an `ω::Float64` for the temporal frequency of the wave
"

# ╔═╡ 4026eb54-9b01-48ba-8bb2-3b2901a65558
ω₀ = 10.0

# ╔═╡ 6104f60d-d961-4963-b7b9-b8a6f5c891be
k₀ = ω₀ # ω/c (c=1)

# ╔═╡ cb2f3944-438f-4a42-81d7-00f9ddb19a67
ω = 1.1ω₀

# ╔═╡ bcaf5797-af13-443c-b670-84ea72d46f8d
md"
#### Unit cell dimensions and discretization
For each periodic dimension, give a number of grid points per unit cell and a size of 
the unit cell: `dims::NTuple{N, Tuple{Int64,Float64}}`
"

# ╔═╡ 1e27a04a-0024-4dd4-8ee2-2cae0ef23d53
Nmodes = 50

# ╔═╡ 017beeeb-327c-4c6b-b60a-e0d5041149d0
λ = 2*pi/k₀

# ╔═╡ df75e240-a9f8-4543-be3f-bf034b152166
md"
#### Wave medium
Planewaves are the eigenmodes of uniform media. Define the medium properties with
either constructor
```julia
struct UniformMedium{T <: Number}
    ϵ::T
    μ::T
end

Vacuum(x=Bool) = UniformMedium(one(x), one(x))
```

The default value is always vacuum
"

# ╔═╡ 22fd2b7b-4418-4cd6-9ff0-eecdf1775d9b
md"
#### Polarization
This isn't stored in `PlanewaveModes`, but it is necessary information for 2D
photonic crystals with diagonal conductivity matrices, when the TE and TM
polarizations decouple. Tell the solver whether to solve `TE()` or `TM()` when
`N=1`. `Coupled()` exists for `N=2` (to be implemented).
"

# ╔═╡ 324aec31-05f2-4fb2-8825-79b2402ec7ee
pol = TM()

# ╔═╡ ce61314c-5ba5-4422-88ce-21f2b20e8110
md"
#### Incident mode amplitudes
Defining a distribution over the incident planewave modes allows for a solution of the
scattering problem. Two distributions must be specified: one at each port of the device.
"

# ╔═╡ d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
### choose which mode to scatter
modeN = 1

# ╔═╡ 4a5de2d0-28e5-4963-b9b3-95fd13f2dd95
I₁ = [n == modeN ? 1 : 0 for n in 1:Nmodes]

# ╔═╡ 1bc23e37-fb20-4d00-95b4-aa7334e18578
I₂ = zeros(Nmodes)

# ╔═╡ 3e0209fe-f316-436d-b017-422844535579
md"
### Defining the scattering structure
The abstract type `RCWAScatterer{N}` is the main respresentation of a unit cell 
in a periodic scattering structure in `DeltaRCWA`.
The type parameter `N` is an integer to specify that a scatterer has a unit cell
with `N` periodic dimensions, which is the dimensionality of the problem.
Any subtype of this can be used to define a problem, but the useful subtypes of it
are `RCWASheet{N} <: RCWAScatterer{N}` which implements an interface for
concrete scatterers with variable surface impedances, and concrete subtypes of
`RCWAStack{N} <: RCWAScatterer{N}` exported by `DeltaRCWA` which represent
sequences of sheets.
For variable volume impedances, which would require a more computationally-expensive
solver for the eigenmodes of Maxwell's equations in that structure such as
[S⁴](http://www.stanford.edu/group/fan/S4/), a type `RCWASlab{N} <: RCWAScatterer{N}`
is defined and could be extended with methods to calculate their scattering matrices.

#### Creating a `RCWASheet`
To create a sheet as part of a 2D photonic crystal, create a struct that is a subtype
of `RCWASheet{1}` and store all the geometric/metasurface parameters you need in
your struct to define the electric and magnetic impedances in the unit cell for it.
"

# ╔═╡ a3035662-3261-43c7-a6cb-2cae4c8b8b0f
struct TrivialSheet{N} <: RCWASheet{N} end

# ╔═╡ a4e6f529-b914-4445-a6ee-7a316e259b9e
struct ComplexExpSheet{N, T} <: RCWASheet{N}
    θ::T # incidence angle
    θᵗ::T # transmission angle
	k::T # wavenumber
    d::T # 
end

# ╔═╡ 77281f06-fcc9-4f75-a415-217cf9e581c1
θ = -π/2.0

# ╔═╡ 1192e390-1f87-486b-b6b4-7447f07f407b
θᵗ = -π/3.0

# ╔═╡ 83eb9d5f-c836-4dd1-be15-1320d94cd743
d = cos(θᵗ)-cos(θ)

# ╔═╡ cef19064-86de-4e25-a5cc-4da77f4b6c2c
L = λ/abs(d) # unit-cell width/periodicity of ComplexExpSheet

# ╔═╡ 4125d1a4-2a57-431a-b7ea-ab8f44994143
dims = ((Nmodes, L), )

# ╔═╡ 6bfcdf25-1f4e-4724-a4bc-14a0222e2c2d
sheet = ComplexExpSheet{1, Float64}(θ, θᵗ, ω, d)

# ╔═╡ a13facc7-06b6-4559-bf2f-a43afb3ffcfe
md"
#### Creating a `SheetStack`
In order to put 2 or more sheets in a sequence, separated by a uniform medium
(the one defined in `PlanewaveModes`) create an instance of a `SheetStack`.
Create a Tuple of the sheets you want to scatter off of and create a second Tuple
with the size of the Vacuum gap that separates each of the sheets and pass these to
the `SheetStack` constructor. Note that there is one gap fewer than the number of sheets.
```julia
struct SheetStack{N, L, T<:Tuple{RCWASheet{N}, Vararg{RCWASheet{N}, L}}} <: RCWAStack{N}
    sheets::T
    depths::Tuple{Vararg{Float64, L}}
end
```
"

# ╔═╡ 8297d683-f888-419c-812b-294397b16126
nsheets = 2

# ╔═╡ 989a6a21-71a8-4991-8a55-eb7f1d3d5f70
gap(x) = 2L

# ╔═╡ 73d6e43f-214b-47ea-b2fd-df60dc412ae1
stack = SheetStack(
	Tuple(sheet for i in 1:nsheets),
	Tuple(gap(i) for i in 1:(nsheets-1)),
)

# ╔═╡ af78faae-2d92-4035-81d5-9856d662e876
md"
#### Methods that a `RCWASheet` implements
Any `RCWASheet` may implement the following 12 methods
`σₑˣˣ, σₑˣʸ, σₑʸˣ, σₑʸʸ, σₘˣˣ, σₘˣʸ, σₘʸˣ, σₘʸʸ, ρₑˣˣ, ρₑˣʸ, ρₑʸˣ, ρₑʸʸ`
one for each component of the electric and magnetic conductivity/resistivity matrices.
(We do not provide `ρₘ` because `ρₘ=0` due to free magnetic charge is non-physical.)
These methods return the value of that admittance/impedance component at points in the
unit cell, thereby realizing the metasurface geometry from the parameters.
By default/fallback, these methods all return zero
(thus are either perfect insulators, `σ`, or are perfect conductors, `ρ`)
unless you
[extend `DeltaRCWA`'s method](https://docs.julialang.org/en/v1/manual/modules/#using-and-import-with-specific-identifiers,-and-adding-methods)
to dispatch on your type.
In a 1D problem, the components used depend on the polarization of light.
Currently the implementation in 1D assumes the off-diagonal components are zero and
neglects them.

See the examples/documentation/help for the method signature
"

# ╔═╡ 030ef40d-043c-4e23-9fc4-7a5321cf67f9
function DeltaRCWA.ρₑˣˣ(sheet::ComplexExpSheet, x⃗::Tuple)
	-0.5sin(sheet.θ)*(1-exp(1im*sheet.k*sheet.d*x⃗[1]))
end

# ╔═╡ 6046fe03-6c7f-4497-a1e1-1368ab043489
function DeltaRCWA.σₘʸʸ(sheet::ComplexExpSheet, x⃗::Tuple)
	-2sin(sheet.θ)*(1+exp(1im*sheet.k*sheet.d*x⃗[1]))
end

# ╔═╡ f2332c95-e784-41c6-a8b9-96df60bb58a6
md"
#### Impedance style
To choose between `ρₑ` or `σₑ`, you specify an `ImpedanceStyle` for your types:
- `Impedanceρₑσₘ` (the default for all `RCWASheet`s)
- `Impedanceσₑσₘ` (optionally available for `RCWASheet{T, 1}`)
which lets you choose whether you define the components of the conductivity matrix σₑ
or its inverse so that you can solve a problem with either a perfect electric
conductor or a perfect electric insulator.
This is a trait for your type.

Uncomment and run the code below to reproduce the same example with the other
impedance style.
"

# ╔═╡ bf0fad44-d478-48c6-9c05-33cd6bf61a7b
# begin
# 	# DeltaRCWA.ImpedanceStyle(::ComplexExpSheet) = Impedanceσₑσₘ()
# 	function DeltaRCWA.σₑˣˣ(sheet::ComplexExpSheet, x⃗)
# 		# this is 1 ./ ρₑˣˣ except with a small perturbation to be nonsingular
# 		2 / (-sin(sheet.θ)*(1e-10 + 1-exp(1im*sheet.k*sheet.d*x⃗[1])))
# 	end
# end

# ╔═╡ f89cf740-818a-45fc-8ec2-5c747fc65fa9
md"
### Problem and solution
Create `DeltaRCWAProblem`s, `solve` them, and receive a `DeltaRCWASolution`.
"

# ╔═╡ aad04ca7-b150-49c2-8a8a-4dc792959649
prob = DeltaRCWAProblem(sheet, dims, ω, pol, I₁, I₂)

# ╔═╡ 786a6947-7082-4902-a625-8be4bd3e30e7
begin
	### Display the magnetic conductivity / electric resistivity  along sheet
	ρ = [-0.5sin(sheet.θ)*(1-exp(1im*sheet.k*sheet.d*x)) for x in prob.modes.x⃗[1]]
	σ = [-2sin(sheet.θ)*(1+exp(1im*sheet.k*sheet.d*x)) for x in prob.modes.x⃗[1]]
	plot(prob.modes.x⃗[1],  real.(ρ), label="Re(ρₑˣˣ)")
	plot!(prob.modes.x⃗[1], imag.(ρ), label="Im(ρₑˣˣ)", ls=:dash)
	plot!(prob.modes.x⃗[1], real.(σ), label="Re(σₘʸʸ)")
	plot!(prob.modes.x⃗[1], imag.(σ), label="Im(σₘʸʸ)", ls=:dash)
end

# ╔═╡ 386f61fa-ab15-40bc-b9be-593622ad42da
sol = solve(prob, method=smatrixLinearMap)

# ╔═╡ fe8b8e5a-7de9-4655-ab19-076d58ce0143
stackprob = DeltaRCWAProblem(stack, dims, ω, pol, I₁, I₂)

# ╔═╡ 20e74b68-23a4-47f1-acea-0b244efd25ae
stacksol = solve(stackprob, method=smatrix)

# ╔═╡ 2bc3991d-9626-4464-a532-5ab7705e454b
md"
### Analysis
The `DeltaRCWASolution` objects can have recipes to analyze and visualize them.
#### Plotting
`plot(::DeltaRCWASolution{T1, T2, 1} where {T1, T2})` should just work!
"

# ╔═╡ 77400f50-4e25-4fe6-8a0a-16f6cf6cb150
plot(sol, combine=true, aspect_ratio=1)

# ╔═╡ d8d2123c-92ce-422e-bb68-bb09e689d44c
plot(stacksol, combine=true, clim=(-1, 1), aspect_ratio=1)

# ╔═╡ 10442fb5-d595-4ef6-a1e9-4d588609dbaf
md"""
## Example of 3D interface
"""

# ╔═╡ 4e412a25-1176-4db9-955b-130fd87e734b
Mmodes = 1 # additional points along y dimension

# ╔═╡ feb92d68-8034-44bf-aba1-bec640315d8f
NMmodes = Nmodes * Mmodes

# ╔═╡ 38dca2e1-2788-4bd9-ad9f-1cfb78bab06c
I₁ˣ = [n == m == modeN ? 0 : 0 for n in 1:Nmodes, m in 1:Mmodes]

# ╔═╡ b6fed619-11f8-40a2-80dd-123adaf9c77c
I₂ˣ = [n == m == modeN ? 0 : 0 for n in 1:Nmodes, m in 1:Mmodes]

# ╔═╡ 3aafea78-513f-4a5c-8a1d-5b689a358a4e
I₁ʸ = [n == m == modeN ? 1 : 0 for n in 1:Nmodes, m in 1:Mmodes]

# ╔═╡ 213628fd-77a1-4753-bc30-f78c36ddd801
I₂ʸ = [n == m == modeN ? 0 : 0 for n in 1:Nmodes, m in 1:Mmodes]

# ╔═╡ b4c64a24-8d14-4054-b6d6-ad7952425115
sheet3d = ComplexExpSheet{2, Float64}(θ, θᵗ, ω, d)

# ╔═╡ bf01b9d0-4437-422f-a769-05955a41e2b1
dims3d = ((Nmodes, L), (Mmodes, L/L))

# ╔═╡ c68fc7c6-5eae-4984-82bf-4368ad910610
prob3d = DeltaRCWAProblem(sheet3d, dims3d, ω, Coupled(), hcat(I₁ˣ, I₁ʸ), hcat(I₂ˣ, I₂ʸ))

# ╔═╡ b966638a-fca0-4431-ac6a-e094c03b138c
sol3d = solve(prob3d)

# ╔═╡ 7eb3f960-c89a-40c7-9c20-f157a2c394a1
z⃗ = range(-L, L, length=2Nmodes)

# ╔═╡ d487aace-25f6-4706-985a-c52bb5799612
z⃗₁ = z⃗[z⃗ .< 0]

# ╔═╡ 458e253d-f1f0-4e26-bd4a-831a44cd6722
z⃗₂ = z⃗[z⃗ .> 0]

# ╔═╡ d5d5ac2a-c94d-4695-8c19-63b9147c29fe
begin
	# quick plotting
	CO₁ = rotr90(DeltaRCWA.bfft(exp.(-sol.modes.kz * transpose(im * z⃗₁)) .* sol3d.O₁[NMmodes.+(1:Nmodes)], 1))
	CI₂ = rotr90(DeltaRCWA.bfft(exp.(-sol.modes.kz * transpose(im * z⃗₂)) .* sol3d.I₂[NMmodes.+(1:Nmodes)], 1))
	CI₁ = rotr90(DeltaRCWA.bfft(exp.( sol.modes.kz * transpose(im * z⃗₁)) .* sol3d.I₁[NMmodes.+(1:Nmodes)], 1))
	CO₂ = rotr90(DeltaRCWA.bfft(exp.( sol.modes.kz * transpose(im * z⃗₂)) .* sol3d.O₂[NMmodes.+(1:Nmodes)], 1))
	heatmap(sol.modes.x⃗, z⃗, real.(cat(CI₁ + CO₁, CI₂ + CO₂; dims=1)), xguide="x", yguide="z", aspect_ratio=:equal,color=:RdBu,clim=(-1.0,1.0))
end

# ╔═╡ 8bb7b1a2-ab66-41ac-86d5-6a29a06b2edc
stack3d = SheetStack(
	Tuple(sheet3d for i in 1:nsheets),
	Tuple(gap(i) for i in 1:(nsheets-1)),
)

# ╔═╡ 301ad780-83a1-46f3-a027-ef8a4ca198f6
prob3dstack = DeltaRCWAProblem(stack3d, dims3d, ω, Coupled(), hcat(I₁ˣ, I₁ʸ), hcat(I₂ˣ, I₂ʸ))

# ╔═╡ f155c239-2001-4fd1-84b2-558827c7cdd2
sol3dstack = solve(prob3dstack)

# ╔═╡ d7372f04-0ede-42ee-b9f5-0e96a2004d68
begin
	# quick plotting
	DO₁ = rotr90(DeltaRCWA.bfft(exp.(-sol.modes.kz * transpose(im * z⃗₁)) .* sol3dstack.O₁[NMmodes.+(1:Nmodes)], 1))
	DI₂ = rotr90(DeltaRCWA.bfft(exp.(-sol.modes.kz * transpose(im * z⃗₂)) .* sol3dstack.I₂[NMmodes.+(1:Nmodes)], 1))
	DI₁ = rotr90(DeltaRCWA.bfft(exp.( sol.modes.kz * transpose(im * z⃗₁)) .* sol3dstack.I₁[NMmodes.+(1:Nmodes)], 1))
	DO₂ = rotr90(DeltaRCWA.bfft(exp.( sol.modes.kz * transpose(im * z⃗₂)) .* sol3dstack.O₂[NMmodes.+(1:Nmodes)], 1))
	heatmap(sol.modes.x⃗, z⃗, real.(cat(DI₁ + DO₁, DI₂ + DO₂; dims=1)), xguide="x", yguide="z", aspect_ratio=:equal,color=:RdBu,clim=(-1.0,1.0))
end

# ╔═╡ Cell order:
# ╠═628d5d5d-2753-47e3-a02b-21b4a89a159e
# ╠═9c58eb1a-313f-4f89-9958-7f33c7a072d3
# ╠═55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
# ╠═65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
# ╟─0a627414-a5af-4fc3-852f-c98105e4d860
# ╟─f0ca02cb-c7b0-42d8-b598-98499aa09f12
# ╠═4026eb54-9b01-48ba-8bb2-3b2901a65558
# ╠═6104f60d-d961-4963-b7b9-b8a6f5c891be
# ╠═cb2f3944-438f-4a42-81d7-00f9ddb19a67
# ╟─bcaf5797-af13-443c-b670-84ea72d46f8d
# ╠═1e27a04a-0024-4dd4-8ee2-2cae0ef23d53
# ╠═017beeeb-327c-4c6b-b60a-e0d5041149d0
# ╠═cef19064-86de-4e25-a5cc-4da77f4b6c2c
# ╠═4125d1a4-2a57-431a-b7ea-ab8f44994143
# ╟─df75e240-a9f8-4543-be3f-bf034b152166
# ╟─22fd2b7b-4418-4cd6-9ff0-eecdf1775d9b
# ╠═324aec31-05f2-4fb2-8825-79b2402ec7ee
# ╟─ce61314c-5ba5-4422-88ce-21f2b20e8110
# ╠═d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
# ╠═4a5de2d0-28e5-4963-b9b3-95fd13f2dd95
# ╠═1bc23e37-fb20-4d00-95b4-aa7334e18578
# ╟─3e0209fe-f316-436d-b017-422844535579
# ╠═a3035662-3261-43c7-a6cb-2cae4c8b8b0f
# ╠═a4e6f529-b914-4445-a6ee-7a316e259b9e
# ╠═77281f06-fcc9-4f75-a415-217cf9e581c1
# ╠═1192e390-1f87-486b-b6b4-7447f07f407b
# ╠═83eb9d5f-c836-4dd1-be15-1320d94cd743
# ╠═6bfcdf25-1f4e-4724-a4bc-14a0222e2c2d
# ╟─a13facc7-06b6-4559-bf2f-a43afb3ffcfe
# ╠═8297d683-f888-419c-812b-294397b16126
# ╠═989a6a21-71a8-4991-8a55-eb7f1d3d5f70
# ╠═73d6e43f-214b-47ea-b2fd-df60dc412ae1
# ╟─af78faae-2d92-4035-81d5-9856d662e876
# ╠═030ef40d-043c-4e23-9fc4-7a5321cf67f9
# ╠═6046fe03-6c7f-4497-a1e1-1368ab043489
# ╠═786a6947-7082-4902-a625-8be4bd3e30e7
# ╟─f2332c95-e784-41c6-a8b9-96df60bb58a6
# ╠═bf0fad44-d478-48c6-9c05-33cd6bf61a7b
# ╟─f89cf740-818a-45fc-8ec2-5c747fc65fa9
# ╠═aad04ca7-b150-49c2-8a8a-4dc792959649
# ╠═386f61fa-ab15-40bc-b9be-593622ad42da
# ╠═fe8b8e5a-7de9-4655-ab19-076d58ce0143
# ╠═20e74b68-23a4-47f1-acea-0b244efd25ae
# ╟─2bc3991d-9626-4464-a532-5ab7705e454b
# ╠═77400f50-4e25-4fe6-8a0a-16f6cf6cb150
# ╠═d8d2123c-92ce-422e-bb68-bb09e689d44c
# ╟─10442fb5-d595-4ef6-a1e9-4d588609dbaf
# ╠═4e412a25-1176-4db9-955b-130fd87e734b
# ╠═feb92d68-8034-44bf-aba1-bec640315d8f
# ╠═38dca2e1-2788-4bd9-ad9f-1cfb78bab06c
# ╠═b6fed619-11f8-40a2-80dd-123adaf9c77c
# ╠═3aafea78-513f-4a5c-8a1d-5b689a358a4e
# ╠═213628fd-77a1-4753-bc30-f78c36ddd801
# ╠═b4c64a24-8d14-4054-b6d6-ad7952425115
# ╠═bf01b9d0-4437-422f-a769-05955a41e2b1
# ╠═c68fc7c6-5eae-4984-82bf-4368ad910610
# ╠═b966638a-fca0-4431-ac6a-e094c03b138c
# ╠═7eb3f960-c89a-40c7-9c20-f157a2c394a1
# ╠═d487aace-25f6-4706-985a-c52bb5799612
# ╠═458e253d-f1f0-4e26-bd4a-831a44cd6722
# ╠═d5d5ac2a-c94d-4695-8c19-63b9147c29fe
# ╠═8bb7b1a2-ab66-41ac-86d5-6a29a06b2edc
# ╠═301ad780-83a1-46f3-a027-ef8a4ca198f6
# ╠═f155c239-2001-4fd1-84b2-558827c7cdd2
# ╠═d7372f04-0ede-42ee-b9f5-0e96a2004d68
