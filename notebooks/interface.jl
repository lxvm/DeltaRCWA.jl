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
using LinearAlgebra: Diagonal

# ╔═╡ adbd86c3-f970-4681-bcde-ddda1050eefd
using FFTW: fftfreq, fft, ifft

# ╔═╡ 67fb8117-7d5b-4536-9e36-7dda36997dff
using BlockArrays: Block

# ╔═╡ 65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
using Plots

# ╔═╡ 9bf1be41-50e7-4efe-a047-66ab1a59bc9a
md"
# Replicating tests

Given what Luke has I will create a scatterer that is similar
(except for the difference in convention regarding conductivities vs impedances).
I will also give this material these properties for the TM polarization,
since this would cause infinite electric conductivity but finite magnetic.
"

# ╔═╡ e92ad8d4-9d2d-44cf-8c6e-37f3d6f52c75
# numerically verify the inversion with the Fourier basis transformation
# fft(ifft(inv(ifft(fft(Matrix(Diagonal(N₀(modes.x⃗...)/2)), 2), 1)),2),1) ≈ Matrix(Diagonal(σₑˣˣ(sheet, modes.x⃗)))

# ╔═╡ 060ff5a3-d8d4-4c4d-8295-4d9fcd07b406
# inv(Matrix(Diagonal(M₀(modes.x⃗...)/2))) ≈ Matrix(Diagonal(inv.(M₀(modes.x⃗...)/2)))

# ╔═╡ 5168b508-df57-40bd-b001-32e2d891d6ae
# plot(inv(Matrix(Diagonal(M₀(modes.x⃗...)/2))))

# ╔═╡ 675e8e41-0aac-404e-9ccf-dfe26710fcfd
# plot(Matrix(Diagonal(σₑˣˣ(sheet, modes.x⃗))))

# ╔═╡ 8321e7a9-f8bc-4334-be6a-19937c026b4a
# methods(smatrix)

# ╔═╡ 4a030218-3def-4122-8433-94b6adb95243
# @which smatrix(sheet, modes, pol)

# ╔═╡ a260a293-bcb4-483c-ba20-2af3e3cbfa58
# plot(real.(sol.I₁))
# plot(real.(Luke_sol.I₁))
# plot(real.(sol.O₂))
# plot(real.(ifft(sol.O₂)))
# plot(real.(Luke_sol.O₂))
# plot(real.(ifft(Luke_sol.O₂)))
# plot(real.(Luke_sol.O₁))
# plot(real.(ifft(Luke_sol.O₁)))

# ╔═╡ d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
modeN = 0

# ╔═╡ aa87c7bd-4d9e-49da-8cf0-4a0fb6081588
begin
#### Luke's solver
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
end

# ╔═╡ c34b5661-e650-4fc1-9dcf-285c39ddb983
plt1

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

function σₑˣˣ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
    2 ./ N₀(x⃗...)
    # ComplexF64[1 - im*sin(sheet.k*sheet.d*x[1])/(1+cos(sheet.k*sheet.d*x[1])) for x in Iterators.product(x⃗...)]
end

function σₘʸʸ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
    2M₀(x⃗...)
end

ω = k
sheet = ComplexExpSheet(θ, θᵗ, k, d, L)
modes = PlanewaveModes(ω, ((length(nvec), sheet.L), ), Vacuum())
pol = TM()
prob = DeltaRCWAProblem(sheet, modes, pol, u_p, u_n)
sol = solve(prob)
end

# ╔═╡ a1174ec5-c1e0-4f18-b93e-1d4b93db7a25
β ≈ modes.kz

# ╔═╡ 65cb38cd-89bc-4254-9ebd-c4d583d8cfb7
plot(sol; method=:fft)

# ╔═╡ 786a6947-7082-4902-a625-8be4bd3e30e7
begin
	### Compare the resitivities/conductivities
	plot(xvec,  real.(M₀(modes.x⃗...)), label="Re(M₀)")
	plot!(xvec, imag.(M₀(modes.x⃗...)), label="Im(M₀)")
	# plot!(xvec, real.(2 ./ M₀(sheet, mode_basis.x⃗)), label="Re(2/M₀)")
	# plot!(xvec, imag.(2 ./ M₀(sheet, mode_basis.x⃗)), label="Im(2/M₀)", ylims=(-15, 15))
	plot!(xvec, real.(σₑˣˣ(sheet, modes.x⃗)), label="Re(σₑˣˣ)")
	plot!(xvec, imag.(σₑˣˣ(sheet, modes.x⃗)), label="Im(σₑˣˣ)", ylims=(-5, 5))
end

# ╔═╡ 3f69f8c5-4d79-47f4-973e-3286cc9b4f6a
begin
	plot(xvec, real.(N₀(modes.x⃗...)), label="Re(N₀)")
	plot!(xvec, imag.(N₀(modes.x⃗...)), label="Im(N₀)")
	plot!(xvec, real.(σₘʸʸ(sheet, modes.x⃗)), label="Re(σₘʸʸ)")
	plot!(xvec, imag.(σₘʸʸ(sheet, modes.x⃗)), label="Im(σₘʸʸ)")
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


# ╔═╡ 2d1f452b-4d01-4fca-ae65-a864c4afa842
Luke_sol = DeltaRCWASolution(modes, pol, u_p, u_n, u_out_p, u_out_n)

# ╔═╡ 07e9c441-8486-426c-abd3-75e48c76a203
β[argmax(real.(Luke_sol.O₂))]

# ╔═╡ a10c8dee-6db7-475e-89f2-bc3c064522f9
plot(Luke_sol; method=:fft)

# ╔═╡ 14457c14-f80e-434e-a118-1ea3486fb3e5
# plot(fft(1:10))
begin
plop = [i+4j for i = 1:10, j = 1:10]
heatmap(1:10, 11:20, plop)
end

# ╔═╡ a9b12d46-c30d-4fca-9ce5-d3dd7498cf5d
# check that the S-matrix matches
S ≈ smatrix(sheet, modes, pol)

# ╔═╡ ca68a6e9-6697-43c3-a51c-966c17b101cc
# plot the spectrum of the matrix
plot(S)

# ╔═╡ d676e456-9984-4d24-8c8f-fece0edb1ee8
plot(Matrix(smatrix(sheet, modes, pol)))

# ╔═╡ 34bb55a9-3e30-4fe8-b60a-98d84c749dce
σˣˣ = ifft(fft(Matrix(Diagonal(reshape(σₑˣˣ(sheet, modes.x⃗), :))), 2), 1)

# ╔═╡ Cell order:
# ╠═628d5d5d-2753-47e3-a02b-21b4a89a159e
# ╠═9c58eb1a-313f-4f89-9958-7f33c7a072d3
# ╠═b4a594f1-a38c-4ede-9860-d4a5afae15c5
# ╠═55bdcbfc-3e78-4db7-ab05-7a9abe9fd253
# ╠═2a30d4a7-ef57-478d-93ac-b5756b6f3909
# ╠═adbd86c3-f970-4681-bcde-ddda1050eefd
# ╠═67fb8117-7d5b-4536-9e36-7dda36997dff
# ╠═9bf1be41-50e7-4efe-a047-66ab1a59bc9a
# ╠═aa87c7bd-4d9e-49da-8cf0-4a0fb6081588
# ╠═c34b5661-e650-4fc1-9dcf-285c39ddb983
# ╠═a1174ec5-c1e0-4f18-b93e-1d4b93db7a25
# ╠═65231a8c-f54d-4ffc-bf59-9dc4cc33f61a
# ╠═65324c70-07b4-46b8-9d6f-3b7fc58d3fbf
# ╠═786a6947-7082-4902-a625-8be4bd3e30e7
# ╠═e92ad8d4-9d2d-44cf-8c6e-37f3d6f52c75
# ╠═060ff5a3-d8d4-4c4d-8295-4d9fcd07b406
# ╠═5168b508-df57-40bd-b001-32e2d891d6ae
# ╠═675e8e41-0aac-404e-9ccf-dfe26710fcfd
# ╠═3f69f8c5-4d79-47f4-973e-3286cc9b4f6a
# ╠═8321e7a9-f8bc-4334-be6a-19937c026b4a
# ╠═4a030218-3def-4122-8433-94b6adb95243
# ╠═65cb38cd-89bc-4254-9ebd-c4d583d8cfb7
# ╠═a260a293-bcb4-483c-ba20-2af3e3cbfa58
# ╠═07e9c441-8486-426c-abd3-75e48c76a203
# ╠═d0d638f3-dc93-48fa-b95d-9fc8b20e22f7
# ╠═3eb5ee6e-8581-41d5-9b72-33d82bad4c5b
# ╠═fc55344b-f0ff-44da-812c-43255aa4e7fd
# ╠═51f5b8a1-f5ec-4b42-a2b3-c3bb9f33e063
# ╠═36afdf48-e275-451b-a095-f2e5e57e5167
# ╠═903b1a7b-8f92-451f-92d7-7e8a72a8c7aa
# ╠═2d1f452b-4d01-4fca-ae65-a864c4afa842
# ╠═a10c8dee-6db7-475e-89f2-bc3c064522f9
# ╠═14457c14-f80e-434e-a118-1ea3486fb3e5
# ╠═a9b12d46-c30d-4fca-9ce5-d3dd7498cf5d
# ╠═ca68a6e9-6697-43c3-a51c-966c17b101cc
# ╠═d676e456-9984-4d24-8c8f-fece0edb1ee8
# ╠═34bb55a9-3e30-4fe8-b60a-98d84c749dce
