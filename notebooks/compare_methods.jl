### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 5c283082-518b-45c6-907a-e80945b1a4e4
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ e8668469-d378-463e-a8be-324f2362d90b
using FFTW

# ╔═╡ 189e4b15-a8bf-4b9b-8832-1e41be5ade12
using LinearAlgebra

# ╔═╡ e0b6e55a-d6f3-4125-9796-a8463c7e00cd
using Plots

# ╔═╡ 8f32367b-b19b-4c18-80ff-5edaabcfb0d9
include("NystromMethodQP.jl");

# ╔═╡ a8253062-df3c-11eb-26e1-0dcff18bf886
md"
# Method comparison: single-sheet scattering

This notebook will compare an RCWA-based method with a boundary-integral equation
method.

For references on the metasurface impedance/conductivity parameters see
- [Pérez-Arancibia et al., 2018, Sideways adiabaticity: beyond ray optics for slowly varying metasurfaces](https://doi.org/10.1364/OE.26.030202)
- [Yu et al., 2011, Light Propagation with Phase Discontinuities: Generalized Laws of Reflection and Refraction](https://doi.org/10.1126/science.1210713)
"

# ╔═╡ 7ead1c90-772b-40fe-8388-9f6347316c2a
md"
## Impedance parameters
"

# ╔═╡ 4d714c5c-cbde-4750-83c2-11ea89bf3041
θⁱ = π/2 # incidence angle, measured counter-clockwise from +x direction

# ╔═╡ 460673fc-b9be-4779-b7c9-5f158e54f0c4
θᵗ = π/3 # desired transmitted field angle

# ╔═╡ dcee0976-8f6d-4dba-9d6c-a91561870419
d = cos(θᵗ)-cos(θⁱ)

# ╔═╡ f09b7267-b5a6-46dd-9ff7-ac48f61531b2
k₀ = 10.0 # wavenumber that attains the metasurface design goals

# ╔═╡ 0ac33952-951c-45be-beb3-e0a43f8b6d60
λ₀ = 2π/k₀ # wavelength that attains the metasurface design goals

# ╔═╡ 2eb8bc93-f6e8-4ec8-a013-dad617a13359
L₀ = λ₀/abs(d) # period of metasurface

# ╔═╡ ec08bb0a-cee7-47e6-a6eb-440673799d7c
# Defining susceptibilities
η(x) = sin(abs(θⁱ)).*(1+exp.(im*k₀*d*x))
# η(x) = 0

# ╔═╡ 18f474f7-72be-4a36-afb4-147592bd68df
μ(x) = sin(abs(θⁱ)).*(1-exp.(im*k₀*d*x))
# μ(x) = 0

# ╔═╡ fd018c1d-6a0d-42bb-8605-a599c70709ee
md"
## Incident wave parameters
"

# ╔═╡ 2b2defdf-66a3-44a2-9e78-25738dfa9f02
k = 1.1*k₀ # = ω/c

# ╔═╡ 2a880954-0413-491e-bd55-0b0a9de2e2e4
λ = 2π/k

# ╔═╡ 2af83733-4752-467d-93cc-1db8209d71dd
α = k * cos(θⁱ)

# ╔═╡ 31098999-eba6-474c-87bf-3d8929ee5017
β = Complex(k * sin(θⁱ))

# ╔═╡ 9d9c4e74-83a8-4e89-a118-45e67fb50e8a
# check that the chosen angle of incidence can be represented exactly in the discrete planewave basis

# ╔═╡ 9028ae22-d871-4f56-b48f-35d85bb8b187
# θⁱ = angle(α + im*β) # incidence angle (measured with respect to the x axis)

# ╔═╡ 84b55287-2b2e-41ad-86b2-2dcdd42634e6
md"
## Discretization parameters
"

# ╔═╡ c22c549b-f588-4c13-8540-978af37ed666
N = 100 # number of discretization points / modes along unit cell

# ╔═╡ a6928e0a-b892-4416-ac66-55e28d60aa7e
L = L₀  # Unit cell width (integer multiple of period of metasurface)

# ╔═╡ e9cbab99-d1c7-44d5-a44d-10d47819e055
dx = L/N

# ╔═╡ b2cc77c1-a578-4290-adda-da4b29ca3d5f
x = range(0, step=dx, length=N) # x points in unit cell

# ╔═╡ 1710cfd1-89c3-4bf8-b825-e76e7c898d74
begin
	# plot metasurface conductivity (η), impedance (μ) parameters
	paramPlot = plot(x,real.(η.(x)), ls=:solid, label="Re η")
	paramPlot = plot!(x,imag.(η.(x)), ls=:dash, label="Im η") 
	paramPlot = plot!(x,real.(μ.(x)), ls=:solid, label="Re μ") 
	paramPlot = plot!(x,imag.(μ.(x)), ls=:dash, label="Im μ") 
	plot(paramPlot,lw=3,frame=:box, xguide="x", yguide="η, μ")
end

# ╔═╡ f476676e-b0d6-46fa-829c-0e691abddffd
y = range(-L, L, length=2N) # y points above and below unit cell

# ╔═╡ 6ac3be67-7dae-42e4-95fd-b0df707892bb
y₁ = y[y .< 0] # y points in domain below sheet

# ╔═╡ 116f6639-9fe5-4f6d-92c2-a6441f11a6d2
y₂ = y[y .> 0] # y points in domain above sheet

# ╔═╡ 852b84b4-5112-4a1c-bece-a767cee13768
k⃗ = 2*pi*fftfreq(N, 1/dx) # wavenumbers for Fourier decomposition

# ╔═╡ cd2385cc-0d16-4f21-910e-3a677263b324
i = findfirst(k⃗ .- α .< eps()) # index of incident mode in k⃗

# ╔═╡ 93fc5b34-18cc-4a88-9984-a5762305b946
β⃗ = @. sqrt(Complex(k^2 - k⃗^2))

# ╔═╡ c4e7232c-34a3-4302-b193-1eeeb60ab718
@assert i == findfirst(real.(β⃗ .- β) .< eps())

# ╔═╡ 1ef4e805-3dd9-4f39-a697-ddc0b4f04ad8
md"
## Carlos' solver
"

# ╔═╡ ad90ebc0-99f2-4636-bce8-12abf791a69a
uInc(x,y)= exp(im*α*x+im*β*y)  # incident planewave

# ╔═╡ dd28d0f8-dcf5-45b0-80f7-9aca4c2ca828
# Vertices that define the curve Γ 
ver = [
	L/2   0;   #1
    -L/2   0   #2
]
# In this case, Γ is just a line segment from (L/2,0) to (-L/2,0). 

# ╔═╡ 2e77a377-69af-496e-8e06-59fa8135fd29
# This matrix defines the orientation of the curve
        #starting node  #ending node  #domain on the right #domian on the left  
info = [1               2             1                    0]

# ╔═╡ 2dcd0f66-f906-4358-907b-7aa0a9199c01
# The parameter M defines the number of boundary points to used in the discretization of the BIE
M  = Int(round(max(20,20/λ)))

# ╔═╡ da3f4379-7b2d-4efc-a9b2-5dcccd4a7d40
# Γ contains all the information of the discretized curve that is needed for the numerical solution of the problem
Γ  = SkeletonCctor(ver,info,[M]);

# ╔═╡ d088a142-f43a-4347-92ca-fb0d6569a262
Γ.plt

# ╔═╡ 51170e8f-2802-4450-b0f6-795b5ad9cacd
method = "MK" # Selection of the Nystrom method to be used

# ╔═╡ 31799b9e-1da0-4e02-8899-2cf7e6ce9563
NC = 5      # Number of unit cells considered in the windowed sum (A = NC*L)

# ╔═╡ f8bb8da7-20d0-43da-846e-5130e7b38edc
mat = matricesSkeleton(1,Γ,k,L,α,NC,method); # Construction of the Nystrom matrices

# ╔═╡ 9fc48c4a-f3f4-4a1a-9fb3-b82a08f3c394
x1 = Γ.parts[1].x[:,1] # x-coordinate of the discretization points on Γ

# ╔═╡ df9e0973-7c43-445e-86bb-fa7030c650e1
f1(x) = -2*im*k*η(x)*exp(im*α*x)

# ╔═╡ cf4f4bf4-8a47-401a-b569-905c9b3cd9fd
f2(x) = -2*im*β*exp(im*α*x)

# ╔═╡ b08df8b1-de6f-4b71-a650-48e483d91a3e
# Solving the boundary integral equations
begin
A⁺ = -0.5*I+im*k*Diagonal(η.(x1))*mat.S  # System matrix for Ω⁺
b⁺ = f1.(x1) # boundary data for Ω⁺
ψ⁺ = A⁺\b⁺

A⁻ = -0.5*I+im*k*Diagonal(μ.(x1))*mat.S # System matrix for Ω⁻
b⁻ = f2.(x1) # boundary data for Ω⁻
ψ⁻ = A⁻\b⁻

φ⁺ = 0.5*(ψ⁺+ψ⁻)
φ⁻ = 0.5*(ψ⁺-ψ⁻)
end;

# ╔═╡ ec2d9dea-e4d7-4703-8422-f4d9c4aec3fd
begin
ver⁺ = [L 0;-L 0;-L 2L;L 2L]; 
Ω⁺ = isin(x,y,ver⁺);

ver⁻ = [L 0;-L 0;-L -2L;L -2L]; 
Ω⁻ =isin(x,y,ver⁻);
end;

# ╔═╡ e9847c59-f98f-419d-88c1-3200ec3e3b4a
UInc = reshape(uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx);

# ╔═╡ 93c2a72c-629f-4512-a048-8245ed2ebeb0
# solve for upper potential (~30 seconds)
pot⁺ = potentialsSkeleton(Ω⁺.pts,Γ,k,L,α,NC);

# ╔═╡ 96f0ddeb-d722-4379-ba89-535f8ea15844
begin
	U⁺scat = reshape(pot⁺.SL*φ⁺,Ω⁺.Ny,Ω⁺.Nx)
	U⁺ = Ω⁺.In .* (U⁺scat .+ UInc)
end;

# ╔═╡ 25a97227-9531-40d4-b354-32663e016ad5
Ω⁺mask = isfinite.(Ω⁺.In);

# ╔═╡ 85a067b0-24e1-47e6-87e7-6850f981cf04
# extract the fields in the upper domain
if θⁱ < 0
	CI₂ = reshape(UInc[Ω⁺mask], length.((y₂, x)))
	CO₂ = reshape(U⁺[Ω⁺mask] .- UInc[Ω⁺mask], length.((y₂, x)))
else
	CI₂ = reshape(UInc[Ω⁺mask] .- UInc[Ω⁺mask], length.((y₂, x)))
	CO₂ = reshape(U⁺[Ω⁺mask], length.((y₂, x)))
end;

# ╔═╡ 46922475-76df-4734-8762-b72dca549b08
# solve for lower potential (~30 seconds)
pot⁻ = potentialsSkeleton(Ω⁻.pts,Γ,k,L,α,NC);

# ╔═╡ d218485c-875a-4d8c-b89b-c967c0f3f10f
begin
	U⁻scat = reshape(pot⁻.SL*φ⁻,Ω⁻.Ny,Ω⁻.Nx)
	U⁻ = Ω⁻.In .* (U⁻scat .+ UInc)
end;

# ╔═╡ 8b34269f-8fdf-44d4-9fa4-ad2e31aa0917
Ω⁻mask = isfinite.(Ω⁻.In);

# ╔═╡ bdf6e2d5-30df-4bb5-a755-b71f29d99b41
# extract the fields in the lower domain
if θⁱ < 0
	CI₁ = reshape(UInc[Ω⁻mask] .- UInc[Ω⁻mask], length.((y₁, x)))
	CO₁ = reshape(U⁻[Ω⁻mask], length.((y₁, x)))
else
	CI₁ = reshape(UInc[Ω⁻mask], length.((y₁, x)))
	CO₁ = reshape(U⁻[Ω⁻mask] .- UInc[Ω⁻mask], length.((y₁, x)))
end;

# ╔═╡ 6045bf4a-e6ac-4411-bd0e-171ee51f3c97
md"
## Luke's solver
"

# ╔═╡ a82986c7-007e-4b6c-9afa-530dc062c5ce
begin
	β̂ = Diagonal(β⃗)
    η̃ = ifft(fft(Matrix(Diagonal(η.(x))), 2), 1)
    μ̃ = ifft(fft(Matrix(Diagonal(μ.(x))), 2), 1)
    A = [
        -β̂-k*η̃    -β̂-k*η̃;
        -β̂-k*μ̃     β̂+k*μ̃
    ]
    B = [
        -β̂+k*η̃    -β̂+k*η̃;
        -β̂+k*μ̃     β̂-k*μ̃
    ]
	S = A\B
	if θⁱ < 0
		I₁ = zeros(N)
		I₂ = [n==i ? 1 : 0 for n in 1:N]
	else
		I₁ = [n==i ? 1 : 0 for n in 1:N]
		I₂ = zeros(N)
	end
	inc = [I₁; I₂]
	scat = S*inc
	O₁ = scat[1:N]
	O₂ = scat[N+1:2N]
end;

# ╔═╡ 5d443572-9d5e-4495-867f-e908d9b56716
md"
## Comparisons

To compare the frequency domain results of Luke's solver to the fields produced by
Carlos' solver, we will have to map the scattered amplitudes onto the grid of Carlos'
solver, possibly correcting for an overall complex normalization factor
(amplitude and phase)

We can find that normalization factor by setting the film to be trivial and comparing
the incident and scattered wave

I will compare 
- incident fields in upper domain (these are just the incident wave)
- incident fields in lower domain (these are both zero)
- scattered fields in upper domain
- scattered fields in lower domain
"

# ╔═╡ bbe3111b-18b3-4832-b636-621567a19632
part = real;

# ╔═╡ a39dd1aa-868d-4819-877f-ccd4b856cb4c
heatmap(x,y,part.(cat(CI₁ + CO₁, CI₂ + CO₂; dims=1)),color=:RdBu,clim=(-1.0,1.0),aspect_ratio=:equal,frame=:box)

# ╔═╡ f667b91a-c07d-44b1-968c-55b80f11eb47
begin
	# quick plotting
	BO₁ = rotr90(bfft(exp.(-β⃗ * transpose(im * y₁)) .* O₁, 1))
	BI₂ = rotr90(bfft(exp.(-β⃗ * transpose(im * y₂)) .* I₂, 1))
	BI₁ = rotr90(bfft(exp.( β⃗ * transpose(im * y₁)) .* I₁, 1))
	BO₂ = rotr90(bfft(exp.( β⃗ * transpose(im * y₂)) .* O₂, 1))
	fields = cat(BI₁ + BO₁, BI₂ + BO₂; dims=1)
	heatmap(x, y, part.(fields), xguide="x", yguide="y", aspect_ratio=:equal,color=:RdBu,clim=(-1.0,1.0))
end

# ╔═╡ 6fc9edc2-c2aa-4be8-9f74-d069de42acb6
phase = -0dx # this is a fudge factor to reduce the apparent far-field error # try -5dx

# ╔═╡ 18ebd13e-88cd-453c-8690-1efc55b35900
# Incident, upper
heatmap(x, y₂, part.(exp(im*phase) * BI₂ - CI₂), aspect_ratio=1, xguide="x", yguide="z", title="$(part)(error incident field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ d9b1970e-f8e3-42a6-b2db-c6d4967bcec2
histogram(reshape(part.(exp(im*phase) * BI₂ - CI₂), :), normalize=true, title="density $(part)(error incident field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ e19bf203-1eac-457f-bd61-1d9bfff9aa88
# scattered, upper
heatmap(x, y₂, part.(exp(phase*im) * BO₂ - CO₂), aspect_ratio=1, xguide="x", yguide="z", title="$(part)(error scattered field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ 3c80c613-78f1-4136-8f6e-dd31b1a60e8e
histogram(reshape(part.(exp(im*phase) * BO₂ - CO₂), :), normalize=true, title="density $(part)(error scattered field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ 980eed0d-3ba6-48a1-b31f-8a1446c96444
# Incident, lower
heatmap(x, y₁, part.(exp(im*phase) * BI₁ - CI₁), aspect_ratio=1, xguide="x", yguide="z", title="$(part)(error incident field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ 620c1306-295e-4ed6-84a6-75fde9e27e66
histogram(reshape(part.(exp(im*phase) * BI₁ - CI₁), :), normalize=true, title="density $(part)(error incident field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ 6aa9e5de-9daf-4b9a-9078-b909a9fc5e03
# scattered, lower
heatmap(x, y₁, part.(exp(im*phase) * BO₁ - CO₁), aspect_ratio=1, xguide="x", yguide="z", title="$(part)(error scattered field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ 0ef24eb1-e9d8-442c-ab8e-03f56735f3c5
md"
Increase resolution to show convergence
"

# ╔═╡ 70e5c1b4-caf4-4760-a8b0-1a509604582e
histogram(reshape(part.(exp(im*phase) * BO₁ - CO₁), :), normalize=true, title="density $(part)(error scattered field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ")

# ╔═╡ Cell order:
# ╟─a8253062-df3c-11eb-26e1-0dcff18bf886
# ╠═5c283082-518b-45c6-907a-e80945b1a4e4
# ╠═e8668469-d378-463e-a8be-324f2362d90b
# ╠═189e4b15-a8bf-4b9b-8832-1e41be5ade12
# ╠═e0b6e55a-d6f3-4125-9796-a8463c7e00cd
# ╟─7ead1c90-772b-40fe-8388-9f6347316c2a
# ╠═4d714c5c-cbde-4750-83c2-11ea89bf3041
# ╠═460673fc-b9be-4779-b7c9-5f158e54f0c4
# ╠═dcee0976-8f6d-4dba-9d6c-a91561870419
# ╠═f09b7267-b5a6-46dd-9ff7-ac48f61531b2
# ╠═0ac33952-951c-45be-beb3-e0a43f8b6d60
# ╠═2eb8bc93-f6e8-4ec8-a013-dad617a13359
# ╠═ec08bb0a-cee7-47e6-a6eb-440673799d7c
# ╠═18f474f7-72be-4a36-afb4-147592bd68df
# ╠═1710cfd1-89c3-4bf8-b825-e76e7c898d74
# ╟─fd018c1d-6a0d-42bb-8605-a599c70709ee
# ╠═2b2defdf-66a3-44a2-9e78-25738dfa9f02
# ╠═2a880954-0413-491e-bd55-0b0a9de2e2e4
# ╠═2af83733-4752-467d-93cc-1db8209d71dd
# ╠═31098999-eba6-474c-87bf-3d8929ee5017
# ╠═9d9c4e74-83a8-4e89-a118-45e67fb50e8a
# ╠═cd2385cc-0d16-4f21-910e-3a677263b324
# ╠═c4e7232c-34a3-4302-b193-1eeeb60ab718
# ╠═9028ae22-d871-4f56-b48f-35d85bb8b187
# ╟─84b55287-2b2e-41ad-86b2-2dcdd42634e6
# ╠═c22c549b-f588-4c13-8540-978af37ed666
# ╠═a6928e0a-b892-4416-ac66-55e28d60aa7e
# ╠═e9cbab99-d1c7-44d5-a44d-10d47819e055
# ╠═b2cc77c1-a578-4290-adda-da4b29ca3d5f
# ╠═f476676e-b0d6-46fa-829c-0e691abddffd
# ╠═6ac3be67-7dae-42e4-95fd-b0df707892bb
# ╠═116f6639-9fe5-4f6d-92c2-a6441f11a6d2
# ╠═852b84b4-5112-4a1c-bece-a767cee13768
# ╠═93fc5b34-18cc-4a88-9984-a5762305b946
# ╟─1ef4e805-3dd9-4f39-a697-ddc0b4f04ad8
# ╠═8f32367b-b19b-4c18-80ff-5edaabcfb0d9
# ╠═ad90ebc0-99f2-4636-bce8-12abf791a69a
# ╠═dd28d0f8-dcf5-45b0-80f7-9aca4c2ca828
# ╠═2e77a377-69af-496e-8e06-59fa8135fd29
# ╠═2dcd0f66-f906-4358-907b-7aa0a9199c01
# ╠═da3f4379-7b2d-4efc-a9b2-5dcccd4a7d40
# ╠═d088a142-f43a-4347-92ca-fb0d6569a262
# ╠═51170e8f-2802-4450-b0f6-795b5ad9cacd
# ╠═31799b9e-1da0-4e02-8899-2cf7e6ce9563
# ╠═f8bb8da7-20d0-43da-846e-5130e7b38edc
# ╠═9fc48c4a-f3f4-4a1a-9fb3-b82a08f3c394
# ╠═df9e0973-7c43-445e-86bb-fa7030c650e1
# ╠═cf4f4bf4-8a47-401a-b569-905c9b3cd9fd
# ╠═b08df8b1-de6f-4b71-a650-48e483d91a3e
# ╠═ec2d9dea-e4d7-4703-8422-f4d9c4aec3fd
# ╠═e9847c59-f98f-419d-88c1-3200ec3e3b4a
# ╠═93c2a72c-629f-4512-a048-8245ed2ebeb0
# ╠═96f0ddeb-d722-4379-ba89-535f8ea15844
# ╠═25a97227-9531-40d4-b354-32663e016ad5
# ╠═85a067b0-24e1-47e6-87e7-6850f981cf04
# ╠═46922475-76df-4734-8762-b72dca549b08
# ╠═d218485c-875a-4d8c-b89b-c967c0f3f10f
# ╠═8b34269f-8fdf-44d4-9fa4-ad2e31aa0917
# ╠═bdf6e2d5-30df-4bb5-a755-b71f29d99b41
# ╠═a39dd1aa-868d-4819-877f-ccd4b856cb4c
# ╟─6045bf4a-e6ac-4411-bd0e-171ee51f3c97
# ╠═a82986c7-007e-4b6c-9afa-530dc062c5ce
# ╠═f667b91a-c07d-44b1-968c-55b80f11eb47
# ╟─5d443572-9d5e-4495-867f-e908d9b56716
# ╠═bbe3111b-18b3-4832-b636-621567a19632
# ╠═6fc9edc2-c2aa-4be8-9f74-d069de42acb6
# ╠═18ebd13e-88cd-453c-8690-1efc55b35900
# ╠═d9b1970e-f8e3-42a6-b2db-c6d4967bcec2
# ╠═e19bf203-1eac-457f-bd61-1d9bfff9aa88
# ╠═3c80c613-78f1-4136-8f6e-dd31b1a60e8e
# ╠═980eed0d-3ba6-48a1-b31f-8a1446c96444
# ╠═620c1306-295e-4ed6-84a6-75fde9e27e66
# ╠═6aa9e5de-9daf-4b9a-9078-b909a9fc5e03
# ╠═0ef24eb1-e9d8-442c-ab8e-03f56735f3c5
# ╠═70e5c1b4-caf4-4760-a8b0-1a509604582e
