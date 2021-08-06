### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 628d5d5d-2753-47e3-a02b-21b4a89a159e
import Pkg

# ╔═╡ 9c58eb1a-313f-4f89-9958-7f33c7a072d3
Pkg.activate(".")

# ╔═╡ 55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
using DeltaRCWA

# ╔═╡ 67fb8117-7d5b-4536-9e36-7dda36997dff
using BlockArrays: Block

# ╔═╡ 65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
using Plots

# ╔═╡ fb39de4e-b436-4dd8-a1b4-0301d75f1e3c
struct NonconductingSheet <: RCWASheet{1} end

# ╔═╡ 414f2e99-226e-41f6-ba05-44539584f71f
smatrix(NonconductingSheet(), PlanewaveModes(12.14, ((10, 1.),), Vacuum()), TE())[Block(2,1)]

# ╔═╡ 5fb3118e-c89f-4e5e-951b-4d7aacda21c7
Nmodes=10

# ╔═╡ 9ac0367a-5924-40fe-a3fd-c8758be6dac2
prob = DeltaRCWAProblem(
	SheetStack((NonconductingSheet(), NonconductingSheet()), (1, )),
	PlanewaveModes(12.14, ((Nmodes, 1.),), Vacuum()),
	TE(),
	ComplexF64[n == 1 ? 1 : 0 for n in 1:Nmodes],
	ComplexF64[n == 1 ? 1 : 0 for n in 1:Nmodes],
	# zeros(ComplexF64, Nmodes),
	# zeros(ComplexF64, Nmodes),
)

# ╔═╡ 34fec8cc-236e-4c97-807d-9ce13502f18e
sol = solve(prob)

# ╔═╡ 197e1720-d390-4816-950e-5f06c66669d4
smatrix(
	SheetStack((NonconductingSheet(), NonconductingSheet()), (1, )),
	PlanewaveModes(12.14, ((Nmodes, 1.),), Vacuum()),
	TE()
)[Block(2,1)]

# ╔═╡ f88832a3-ae32-41bd-8829-8cb34709f970
sol.O₂

# ╔═╡ 9bf1be41-50e7-4efe-a047-66ab1a59bc9a
md"
# Replicating tests

Luke has

```julia
#### Matching Carlos's solver
k = 10.0;
λ = 2*pi/k;
θ = -π/2.0; # incidence angle (measured with respect to the x axis)
α = k*cos(θ);
β = k*sin(θ);
uInc(x,y)= @. exp(1im*α*x+1im*β*y);  # incident planewave

θᵗ = -π/8;    # transmitted field angle
d  = cos(θᵗ)-cos(θ); 
L = 2*(2*π)/(k*abs(d));  # Unit cell width
M₀(x) = @. -sin(θ)*(1+exp(1im*k*d*x));
N₀(x) = @. -sin(θ)*(1-exp(1im*k*d*x));

nvec = 0:99;
dx = L/length(nvec)
xvec = [n*dx-L/2 for n in nvec]
kₓ = 2*pi*fftfreq(length(nvec), 1/dx)
β = @. sqrt(Complex(k^2 - kₓ^2))

u_p = fft(uInc(xvec, 0))/length(nvec)
u_n = zeros(length(nvec))
M = Matrix(Diagonal(M₀(xvec)))
N = Matrix(Diagonal(N₀(xvec)))

A = [-Diagonal(β)-k*ifft(fft(M, 2), 1)    -Diagonal(β)-k*ifft(fft(M, 2), 1);
    -Diagonal(β)-k*ifft(fft(N, 2), 1)     Diagonal(β)+k*ifft(fft(N, 2), 1)]
B = [-Diagonal(β)+k*ifft(fft(M, 2), 1)    -Diagonal(β)+k*ifft(fft(M, 2), 1);
    -Diagonal(β)+k*ifft(fft(N, 2), 1)     Diagonal(β)-k*ifft(fft(N, 2), 1)]
S = A\B

u_in = [u_p; u_n]
u_out = S*u_in
```

I will create a scatterer that is similar (except for the difference in convention
regarding conductivities vs impedances).
I will also give this material these properties for the TM polarization,
since this would cause infinite electric conductivity but finite magnetic.
"

# ╔═╡ 248713b5-fef5-4109-a3e4-d84eec4c071e
begin
	"parametrization of Luke's sheet profile"
	struct ComplexExpSheet{T <: Real} <: RCWASheet{1}
		θ::T
		θᵗ::T
		k::T
		d::T
		L::T
	end

	"Outer Constructor"
	function ComplexExpSheet(θ::T, θᵗ::T, k::T) where T
		d = cos(θᵗ) - cos(θ)
		ComplexExpSheet{T}(θ, θᵗ, k, d, 2*(2*π)/(k*abs(d)))
	end
end

# ╔═╡ ad61805c-93ad-4856-b661-b08d3fe83532
function M₀(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
	ComplexF64[-sin(sheet.θ)*(1+exp(im*sheet.k*sheet.d*x[1])) for x in Iterators.product(x⃗...)]
end

# ╔═╡ e77f6808-1138-47eb-b1ef-f6bef2c58c8f
function N₀(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
	ComplexF64[-sin(sheet.θ)*(1-exp(im*sheet.k*sheet.d*x[1])) for x in Iterators.product(x⃗...)]
end

# ╔═╡ cbc351a6-26fa-40c1-b3a3-d5ab2c4e0250
function σₑˣˣ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
	2 ./ M₀(sheet, x⃗)
	# ComplexF64[1 - im*sin(sheet.k*sheet.d*x[1])/(1+cos(sheet.k*sheet.d*x[1])) for x in Iterators.product(x⃗...)]
end

# ╔═╡ 6fff38a0-17f7-4a2a-8563-d3218f25407f
function σₘʸʸ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
	2N₀(sheet, x⃗)
end

# ╔═╡ 0d4c6c3f-4258-4e73-9df7-929c87b7d06a
1 / 5

# ╔═╡ 928cce69-2599-42eb-a144-e74aec87bf3d
begin
	θ = -π/2.0
	θᵗ= -π/8.0
	d = cos(θᵗ) - cos(θ)
	k =10.0
	ω = k
	L = 2*(2*π)/(k*abs(d))
	α = k*cos(θ)
	β = k*sin(θ)
	sheet = ComplexExpSheet(θ, θᵗ, k)
	nmodes = 100
	I₁ = ComplexF64[n == 1 ? 1 : 0 for n in 1:nmodes]
	I₂ = ComplexF64[0 for n in 1:nmodes]
	mode_basis = PlanewaveModes(ω, ((nmodes, sheet.L), ), Vacuum())
	pol = TM()
end

# ╔═╡ 3373913a-fd89-4379-8a7a-d0a3910d6135
σₑˣˣ(sheet, mode_basis.x⃗)

# ╔═╡ 9b7306d3-4087-42a9-9568-66503454dbce
M₀(sheet, mode_basis.x⃗)

# ╔═╡ e92ad8d4-9d2d-44cf-8c6e-37f3d6f52c75
fft(ifft(inv(ifft(fft(Matrix(Diagonal(M₀(sheet, mode_basis.x⃗)/2)), 2), 1)),2),1)

# ╔═╡ 33550921-61fc-42f5-83ce-f16bea5a7bda
begin
	uInc(x,y)= @. exp(1im*α*x+1im*β*y)
	nvec = 0:99
	dx = L/length(nvec)
	xvec = [n*dx-L/2 for n in nvec]
	u_p = fft(uInc(xvec, 0))/length(nvec)
end

# ╔═╡ 786a6947-7082-4902-a625-8be4bd3e30e7
begin
	plot(xvec, real.(M₀(sheet, mode_basis.x⃗)), label="Re(M₀)")
	plot!(xvec, imag.(M₀(sheet, mode_basis.x⃗)), label="Im(M₀)")
	# plot!(xvec, real.(2 ./ M₀(sheet, mode_basis.x⃗)), label="Re(2/M₀)")
	# plot!(xvec, imag.(2 ./ M₀(sheet, mode_basis.x⃗)), label="Im(2/M₀)", ylims=(-15, 15))
	plot!(xvec, real.(σₑˣˣ(sheet, mode_basis.x⃗)), label="Re(σₑˣˣ)")
	plot!(xvec, imag.(σₑˣˣ(sheet, mode_basis.x⃗)), label="Im(σₑˣˣ)", ylims=(-5, 5))
end

# ╔═╡ 3f69f8c5-4d79-47f4-973e-3286cc9b4f6a
begin
	plot(xvec, real.(N₀(sheet, mode_basis.x⃗)), label="Re(N₀)")
	plot!(xvec, imag.(N₀(sheet, mode_basis.x⃗)), label="Im(N₀)")
	plot!(xvec, real.(σₘʸʸ(sheet, mode_basis.x⃗)), label="Re(σₘʸʸ)")
	plot!(xvec, imag.(σₘʸʸ(sheet, mode_basis.x⃗)), label="Im(σₘʸʸ)")
end

# ╔═╡ a0a6a85d-78d1-40f5-96d9-80241c22e735
test_prob = DeltaRCWAProblem(sheet,	mode_basis, pol, I₁, I₂)

# ╔═╡ 00c1c72e-426f-46c5-a37a-411697d71189
test_sol = solve(test_prob)

# ╔═╡ 2900e92b-4857-4d40-a2dc-9e8bfb61bf72
typeof(sheet)

# ╔═╡ fbb775af-15cc-4d71-af64-75c68e9c8461
σₑˣˣ(sheet, mode_basis.x⃗)

# ╔═╡ 8321e7a9-f8bc-4334-be6a-19937c026b4a
methods(smatrix)

# ╔═╡ 4a030218-3def-4122-8433-94b6adb95243
# @which
smatrix(sheet, mode_basis, pol)

# ╔═╡ 65cb38cd-89bc-4254-9ebd-c4d583d8cfb7
begin
	heatmap(mode_basis.x⃗..., map(x -> x, mode_basis.x⃗...), real.([0 for e in Iterators.product(test_sol.I₁, mode_basis.x⃗...)]))
end

# ╔═╡ 40ea61a7-c43c-45e2-a918-c556172b693f
plot(real.(ifft(exp.(im .* test_sol.modes.kz) .* reshape(ComplexF64[n == 4 ? 1 : 0 for n in 1:nmodes], length.(test_sol.modes.k⃗)))))

# ╔═╡ Cell order:
# ╠═628d5d5d-2753-47e3-a02b-21b4a89a159e
# ╠═9c58eb1a-313f-4f89-9958-7f33c7a072d3
# ╠═55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
# ╠═67fb8117-7d5b-4536-9e36-7dda36997dff
# ╠═fb39de4e-b436-4dd8-a1b4-0301d75f1e3c
# ╠═414f2e99-226e-41f6-ba05-44539584f71f
# ╠═5fb3118e-c89f-4e5e-951b-4d7aacda21c7
# ╠═9ac0367a-5924-40fe-a3fd-c8758be6dac2
# ╠═34fec8cc-236e-4c97-807d-9ce13502f18e
# ╠═197e1720-d390-4816-950e-5f06c66669d4
# ╠═f88832a3-ae32-41bd-8829-8cb34709f970
# ╟─9bf1be41-50e7-4efe-a047-66ab1a59bc9a
# ╠═248713b5-fef5-4109-a3e4-d84eec4c071e
# ╠═ad61805c-93ad-4856-b661-b08d3fe83532
# ╠═e77f6808-1138-47eb-b1ef-f6bef2c58c8f
# ╠═cbc351a6-26fa-40c1-b3a3-d5ab2c4e0250
# ╠═6fff38a0-17f7-4a2a-8563-d3218f25407f
# ╠═3373913a-fd89-4379-8a7a-d0a3910d6135
# ╠═0d4c6c3f-4258-4e73-9df7-929c87b7d06a
# ╠═9b7306d3-4087-42a9-9568-66503454dbce
# ╠═65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
# ╠═786a6947-7082-4902-a625-8be4bd3e30e7
# ╠═e92ad8d4-9d2d-44cf-8c6e-37f3d6f52c75
# ╠═3f69f8c5-4d79-47f4-973e-3286cc9b4f6a
# ╠═33550921-61fc-42f5-83ce-f16bea5a7bda
# ╠═928cce69-2599-42eb-a144-e74aec87bf3d
# ╠═a0a6a85d-78d1-40f5-96d9-80241c22e735
# ╠═00c1c72e-426f-46c5-a37a-411697d71189
# ╠═2900e92b-4857-4d40-a2dc-9e8bfb61bf72
# ╠═fbb775af-15cc-4d71-af64-75c68e9c8461
# ╠═8321e7a9-f8bc-4334-be6a-19937c026b4a
# ╠═4a030218-3def-4122-8433-94b6adb95243
# ╠═65cb38cd-89bc-4254-9ebd-c4d583d8cfb7
# ╠═40ea61a7-c43c-45e2-a918-c556172b693f
