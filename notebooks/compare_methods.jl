### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 5c283082-518b-45c6-907a-e80945b1a4e4
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ f334af2a-897b-42bc-80c6-756e6e3519a8
using Revise

# ╔═╡ 898fa516-591d-4831-8c65-f0e758922d15
using DeltaRCWA

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
# Testing multilayer films

For references on the metasurface impedance/conductivity parameters see
- [Pérez-Arancibia et al., 2018, Sideways adiabaticity: beyond ray optics for slowly varying metasurfaces](https://doi.org/10.1364/OE.26.030202)
- [Yu et al., 2011, Light Propagation with Phase Discontinuities: Generalized Laws of Reflection and Refraction](https://doi.org/10.1126/science.1210713)
"

# ╔═╡ 1ef4e805-3dd9-4f39-a697-ddc0b4f04ad8
md"
## Carlos' solver
"

# ╔═╡ ad90ebc0-99f2-4636-bce8-12abf791a69a
# Defining constants
begin
k = 10.0   # wavenumber k²=ω²μ₀ϵ₀

λ = 2*pi/k # wavelength 

θ = -π/2.0 # incidence angle (measured with respect to the x axis)

α = k*cos(θ)

β = k*sin(θ)

uInc(x,y)= exp(im*α*x+im*β*y)  # incident planewave
end;

# ╔═╡ fea6690d-7602-462b-85dd-4fd504918933
begin
θᵗ = -π/3;    # transmitted field angle

d  = cos(θᵗ)-cos(θ); 

# Defining susceptibilities
η(x) = -sin(θ).*(1+exp.(im*k*d*x));
# η(x) = 0#-sin(θ).*(1+exp.(im*k*d*x));

μ(x) = 1e-8-sin(θ).*(1-exp.(im*k*d*x));
# μ(x) = 1000000000 #-sin(θ).*(1-exp.(im*k*d*x));

L = 2*(2*π)/(k*abs(d));  # Unit cell width
end;

# ╔═╡ da3f4379-7b2d-4efc-a9b2-5dcccd4a7d40
# Defining geometry
begin
# Vertices that define the curve Γ 
ver = [  L/2   0;   #1 
        -L/2   0]  #2 
# In this case, Γ is just a line segment from (L/2,0) to (-L/2,0). 

# This matrix defines the orientation of the curve
        #starting node  #ending node  #domain on the right #domian on the left  
info = [1               2             1                    0]

# The parameter P defines the number of boundary points to used in the discretization of the BIE
P  = Int(round(max(20,20/λ)))

# Γ contains all the information of the discretized curve that is needed for the numerical solution of the problem
Γ  = SkeletonCctor(ver,info,[P])
end;

# ╔═╡ 1710cfd1-89c3-4bf8-b825-e76e7c898d74
begin
	# plot metasurface conductivity (η), impedance (μ) parameters
	xPlot = -L/2:0.001:L/2;
	paramPlot = plot(xPlot,real.(η.(xPlot)), ls=:solid, label="Re η")
	paramPlot = plot!(xPlot,imag.(η.(xPlot)), ls=:dash, label="Re η") 
	paramPlot = plot!(xPlot,real.(μ.(xPlot)), ls=:solid, label="Re μ") 
	paramPlot = plot!(xPlot,imag.(μ.(xPlot)), ls=:dash, label="Re μ") 
	plot(paramPlot,lw=3,frame=:box, xguide="x", yguide="η, μ")
end

# ╔═╡ b08df8b1-de6f-4b71-a650-48e483d91a3e
# Solving the boundary integral equations
begin
method = "MK" # Selection of the Nystrom method to be used

NC  = 5      # Number of unit cells considered in the windowed sum (A = NC*L)

mat = matricesSkeleton(1,Γ,k,L,α,NC,method) # Construction of the Nystrom matrices
    
x1 = Γ.parts[1].x[:,1] # x-coordinate of the discretization points on Γ
    
Id = Matrix(I,Γ.Npts,Γ.Npts) # Indentity matrix

f1(x) = -2*im*k*η(x)*exp(im*α*x)
A⁺ = -0.5*Id+im*k*Diagonal(η.(x1))*mat.S  # System matrix for Ω⁺
b⁺ = f1.(x1) # boundary data for Ω⁺
ψ⁺ = A⁺\b⁺

f2(x) = -2*im*β*exp(im*α*x)
A⁻ = -0.5*Id+im*k*Diagonal(μ.(x1))*mat.S # System matrix for Ω⁻
b⁻ =  f2.(x1) # boundary data for Ω⁻
ψ⁻ = A⁻\b⁻

φ⁺ = 0.5*(ψ⁺+ψ⁻)
φ⁻ = 0.5*(ψ⁺-ψ⁻)
end;

# ╔═╡ c9f9390a-5d26-4791-a023-558e0a911d08
begin
	lims = [-L L -L L]
	h=0.025
	x = lims[1]:h:lims[2]
	y = lims[3]:h:lims[4]
end;

# ╔═╡ ec2d9dea-e4d7-4703-8422-f4d9c4aec3fd
begin
ver⁺ = [L 0;-L 0;-L 100;L 100]; 
Ω⁺ = isin(x,y,ver⁺);

ver⁻ = [L 0;-L 0;-L -100;L -100]; 
Ω⁻ =isin(x,y,ver⁻);
end;

# ╔═╡ f17855ee-a00d-4b39-bf39-4470a25bf32b
begin
	UInc = reshape(uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx)
	@assert UInc ≈ reshape(uInc.(Ω⁻.pts[:,1],Ω⁻.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx)
end

# ╔═╡ 93c2a72c-629f-4512-a048-8245ed2ebeb0
# solve for upper potential (~30 seconds)
pot⁺ = potentialsSkeleton(Ω⁺.pts,Γ,k,L,α,NC);

# ╔═╡ 96f0ddeb-d722-4379-ba89-535f8ea15844
begin
	U⁺scat = reshape(pot⁺.SL*φ⁺,Ω⁺.Ny,Ω⁺.Nx)
	U⁺   =  Ω⁺.In .* (U⁺scat .+ UInc)
end;

# ╔═╡ 46922475-76df-4734-8762-b72dca549b08
# solve for lower potential (~30 seconds)
pot⁻= potentialsSkeleton(Ω⁻.pts,Γ,k,L,α,NC);

# ╔═╡ d218485c-875a-4d8c-b89b-c967c0f3f10f
begin
	U⁻scat = reshape(pot⁻.SL*φ⁻,Ω⁻.Ny,Ω⁻.Nx)
	U⁻   =  Ω⁻.In .* (U⁻scat .+ UInc)
end;

# ╔═╡ a39dd1aa-868d-4819-877f-ccd4b856cb4c
begin
	pltReal =  heatmap(x,y,real.(U⁺),color=:RdBu,clim=(-1.0,1.0))
	pltReal =  heatmap!(x,y,real.(U⁻),color=:RdBu,clim=(-1.0,1.0))
	pltReal = plot!(pltReal,legend=false,aspect_ratio=:equal,frame=:box)
end

# ╔═╡ a022ef8b-03cc-42cc-ac9e-4da96e3edecb
begin
	pltImag =  heatmap(x,y,imag.(U⁺),color=:RdBu,clim=(-1.0,1.0))
	pltImag =  heatmap!(x,y,imag.(U⁻),color=:RdBu,clim=(-1.0,1.0))
	pltImag = plot!(pltImag,legend=false,aspect_ratio=:equal,frame=:box)
end

# ╔═╡ 66a8554e-1b71-4776-a903-c9ee54f1d64b
	heatmap(x,y,abs2.(U⁻),color=:RdBu,clim=(-1.0,1.0))

# ╔═╡ 78e6c13a-ba9a-4b5f-8ae6-51ebeb7f098e
md"
## DeltaRCWA solver
"

# ╔═╡ 6045bf4a-e6ac-4411-bd0e-171ee51f3c97
md"
## Luke's solver
"

# ╔═╡ a82986c7-007e-4b6c-9afa-530dc062c5ce
begin
nvec = 0:99
dx = L/length(nvec)
xvec = [n*dx-L/2 for n in nvec]
kₓ = 2*pi*fftfreq(length(nvec), 1/dx)
βz = @. sqrt(Complex(k^2 - kₓ^2))

u_p = fft(uInc.(xvec, 0))/length(nvec)
u_n = zeros(length(nvec))
M = Matrix(Diagonal(η.(xvec)))
N = Matrix(Diagonal(μ.(xvec)))

A = [-Diagonal(βz)-k*ifft(fft(M, 2), 1)    -Diagonal(βz)-k*ifft(fft(M, 2), 1);
    -Diagonal(βz)-k*ifft(fft(N, 2), 1)     Diagonal(βz)+k*ifft(fft(N, 2), 1)]
B = [-Diagonal(βz)+k*ifft(fft(M, 2), 1)    -Diagonal(βz)+k*ifft(fft(M, 2), 1);
    -Diagonal(βz)+k*ifft(fft(N, 2), 1)     Diagonal(βz)-k*ifft(fft(N, 2), 1)]
S = A\B

u_in = [u_p; u_n]
u_out = S*u_in
u_out_p = u_out[1:length(nvec)]
u_out_n = u_out[length(nvec)+1:2*length(nvec)]
end;

# ╔═╡ 07b0e662-b295-46ff-a744-13b5e6fef561
begin
	struct ComplexExpSheet{T} <: RCWASheet{T, 1}
		θ::T
		θᵗ::T
		k::T
		d::T
		L::T
	end

	### Define how to convert between η/μ and conductivity matrix conventions
	DeltaRCWA.σₑˣˣ(::ComplexExpSheet, x⃗) = 2 ./ μ.(x⃗[1])
	DeltaRCWA.σₘʸʸ(::ComplexExpSheet, x⃗) = 2η.(x⃗[1])

	ω = k
	sheet = ComplexExpSheet(θ, θᵗ, k, d, L)
	dims = ((length(nvec), sheet.L), )
	pol = TM()
	prob = DeltaRCWAProblem(sheet, dims, ω, pol,
		zeros(length(xvec)),
		[n==1 ? 1 : 0 for n in 1:length(nvec)]
		# ifft(uInc.(xvec, 0)), # too noisy to plot
	)
	sol = solve(prob)
end;

# ╔═╡ 3e3b9f60-6791-4922-8d3a-0dbe98437060
plot(sol, method=:fft, combine=false,
	# mask=Bool[n==argmax(abs.(sol.O₂)) ? false : true for n in 1:length(nvec)]
)

# ╔═╡ 5379c325-69fa-4f3e-957d-644d40155997
any(iszero.(prob.modes.kz))

# ╔═╡ 926a554f-e587-40fa-99b5-fad6171fd551
# this is the mode 
prob.modes.k⃗[1][argmax(abs.(sol.O₂))], prob.modes.kz[argmax(abs.(sol.O₂))]

# ╔═╡ 63b2e16e-5614-472d-88d5-09ef17c6aeaa
begin
	# wrangle Luke's solution into DeltaRCWA for plotting
	Lukemodes = PlanewaveModes(k, Vacuum(), ((length(nvec), L), ), prob.modes.x⃗, (kₓ, ), βz)
	Lukesol = DeltaRCWA.DeltaRCWASolution(Lukemodes, pol, u_n, u_p, u_out_n, u_out_p)
	plot(Lukesol)
end

# ╔═╡ d2350cc9-b4a1-4ab2-ab0a-39c0dad37681
md"
## Comparisons

### Luke/DeltaRCWA
Since these are basically the same solver, we can calculate the difference between
these models in terms of the norm of the difference between the scattered amplitudes
"

# ╔═╡ 1fbcf167-e742-48a7-a939-445ac0798550
# error of scattered amplitudes in port 1
norm(u_out_n - sol.O₁)

# ╔═╡ 1b62b57e-d6c8-4fa8-8d4a-f7e62e5b795f
# error of scattered amplitudes in port 2
norm(u_out_p - sol.O₂)

# ╔═╡ 5d443572-9d5e-4495-867f-e908d9b56716
md"
### DeltaRCWA/Carlos

To compare the frequency domain results of DeltaRCWA to the fields produced by
Carlos' solver, we will have the scattered amplitudes onto the grid of Carlos'
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

# ╔═╡ 4a1aa64a-96da-4ecc-b98d-9623da487b9f
begin
	# find the subsets of the x, y points in each domain (Ω⁺, Ω⁻)
	xisinΩ⁺ = [1.0 in Ω⁺.In[:, i] ? true : false for i in eachindex(x)]
	yisinΩ⁺ = [1.0 in Ω⁺.In[i, :] ? true : false for i in eachindex(y)]
	xinΩ⁺ = x[findfirst(xisinΩ⁺):findlast(xisinΩ⁺)]
	yinΩ⁺ = y[findfirst(yisinΩ⁺):findlast(yisinΩ⁺)]
	xisinΩ⁻ = [1.0 in Ω⁻.In[:, i] ? true : false for i in eachindex(x)]
	yisinΩ⁻ = [1.0 in Ω⁻.In[i, :] ? true : false for i in eachindex(y)]
	xinΩ⁻ = x[findfirst(xisinΩ⁻):findlast(xisinΩ⁻)]
	yinΩ⁻ = y[findfirst(yisinΩ⁻):findlast(yisinΩ⁻)]
end;

# ╔═╡ e085b47a-a8a3-4620-9fce-f031c837e33b
begin
	k⃗ = sol.modes.k⃗[1]
	kz = prob.modes.kz
	# Luke's equivalent
	# k⃗ = kₓ
	# kz = βz
end;

# ╔═╡ fb36dabc-4aa0-40dd-ab3b-c7b97f34ab37
begin
	mask = fill(true, size(kz))
	# mask = Bool[n==argmax(abs.(sol.O₂)) ? false : true for n in 1:length(nvec)]
	O₂ = sol.O₂ .* mask
    I₂ = sol.I₂ .* mask
    I₁ = sol.I₁ .* mask
    O₁ = sol.O₁ .* mask
	# To compare with Luke, use these
	# O₂ = Lukesol.O₂ .* mask
	# I₂ = Lukesol.I₂ .* mask
	# I₁ = Lukesol.I₁ .* mask
	# O₁ = Lukesol.O₁ .* mask
end;

# ╔═╡ 039bdb68-c4ce-4a42-ab4a-8f91c5ac29d6
begin
	# DeltaRCWA: compute the fields in the upper domain
	AO₂ = zeros(ComplexF64, length.((xinΩ⁺, yinΩ⁺)))
	AI₂ = zeros(ComplexF64, length.((xinΩ⁺, yinΩ⁺)))
	for j in eachindex(xinΩ⁺)
		for i in eachindex(yinΩ⁺)
			AO₂[j, i] = sum(exp.((yinΩ⁺[i] * im) .*  kz .+ (xinΩ⁺[j] * im) .* k⃗) .* O₂)
			AI₂[j, i] = sum(exp.((yinΩ⁺[i] * im) .* -kz .+ (xinΩ⁺[j] * im) .* k⃗) .* I₂)
		end
	end
	AO₂ = transpose(reverse(AO₂, dims=1))
	AI₂ = transpose(reverse(AI₂, dims=1))
end;

# ╔═╡ 317c15e3-dc43-494a-8ad6-377b1d7a9ff8
begin
	# DeltaRCWA: compute the fields in the lower domain
	AI₁ = zeros(ComplexF64, length.((xinΩ⁻, yinΩ⁻)))
	AO₁ = zeros(ComplexF64, length.((xinΩ⁻, yinΩ⁻)))
	for j in eachindex(xinΩ⁻)
		for i in eachindex(yinΩ⁻)
			AI₁[j, i] = sum(exp.((yinΩ⁻[i] * im) .*  kz .+ (xinΩ⁻[j] * im) .* k⃗) .* I₁)
			AO₁[j, i] = sum(exp.((yinΩ⁻[i] * im) .* -kz .+ (xinΩ⁻[j] * im) .* k⃗) .* O₁)
		end
	end
	AO₁ = transpose(reverse(AO₁, dims=1))
	AI₁ = transpose(reverse(AI₁, dims=1))
end;

# ╔═╡ 85a067b0-24e1-47e6-87e7-6850f981cf04
begin
	# Carlos: extract the fields in the upper domain
	Ω⁺mask = isfinite.(Ω⁺.In)
	CI₂ = reshape(UInc[Ω⁺mask], length.((yinΩ⁺, xinΩ⁺)))
	CO₂ = reshape(U⁺[Ω⁺mask] .- UInc[Ω⁺mask], length.((yinΩ⁺, xinΩ⁺)))
end;

# ╔═╡ bdf6e2d5-30df-4bb5-a755-b71f29d99b41
begin
	# Carlos: extract the fields in the lower domain
	Ω⁻mask = isfinite.(Ω⁻.In)
	CI₁ = reshape(UInc[Ω⁻mask] .- UInc[Ω⁻mask], length.((yinΩ⁻, xinΩ⁻)))
	CO₁ = reshape(U⁻[Ω⁻mask], length.((yinΩ⁻, xinΩ⁻)))
end;

# ╔═╡ bbe3111b-18b3-4832-b636-621567a19632
part = real;

# ╔═╡ 980eed0d-3ba6-48a1-b31f-8a1446c96444
# Incident, upper
heatmap(part.(AI₁ .- CI₁), clim=(-1, 1))

# ╔═╡ e19bf203-1eac-457f-bd61-1d9bfff9aa88
# scattered, upper
heatmap(part.(AO₂ .- CO₂), clim=(-1, 1))

# ╔═╡ 6aa9e5de-9daf-4b9a-9078-b909a9fc5e03
# scattered, lower
heatmap(part.(AO₁ .- CO₁), clim=(-1, 1))

# ╔═╡ 18ebd13e-88cd-453c-8690-1efc55b35900
# Incident, lower
heatmap(part.(AI₂ .- CI₂), clim=(-1, 1))

# ╔═╡ Cell order:
# ╟─a8253062-df3c-11eb-26e1-0dcff18bf886
# ╠═5c283082-518b-45c6-907a-e80945b1a4e4
# ╠═f334af2a-897b-42bc-80c6-756e6e3519a8
# ╠═898fa516-591d-4831-8c65-f0e758922d15
# ╠═e8668469-d378-463e-a8be-324f2362d90b
# ╠═189e4b15-a8bf-4b9b-8832-1e41be5ade12
# ╠═e0b6e55a-d6f3-4125-9796-a8463c7e00cd
# ╟─1ef4e805-3dd9-4f39-a697-ddc0b4f04ad8
# ╠═8f32367b-b19b-4c18-80ff-5edaabcfb0d9
# ╠═ad90ebc0-99f2-4636-bce8-12abf791a69a
# ╠═fea6690d-7602-462b-85dd-4fd504918933
# ╠═da3f4379-7b2d-4efc-a9b2-5dcccd4a7d40
# ╠═1710cfd1-89c3-4bf8-b825-e76e7c898d74
# ╠═b08df8b1-de6f-4b71-a650-48e483d91a3e
# ╠═c9f9390a-5d26-4791-a023-558e0a911d08
# ╠═ec2d9dea-e4d7-4703-8422-f4d9c4aec3fd
# ╠═f17855ee-a00d-4b39-bf39-4470a25bf32b
# ╠═93c2a72c-629f-4512-a048-8245ed2ebeb0
# ╠═96f0ddeb-d722-4379-ba89-535f8ea15844
# ╠═46922475-76df-4734-8762-b72dca549b08
# ╠═d218485c-875a-4d8c-b89b-c967c0f3f10f
# ╠═a39dd1aa-868d-4819-877f-ccd4b856cb4c
# ╠═a022ef8b-03cc-42cc-ac9e-4da96e3edecb
# ╠═66a8554e-1b71-4776-a903-c9ee54f1d64b
# ╟─78e6c13a-ba9a-4b5f-8ae6-51ebeb7f098e
# ╠═07b0e662-b295-46ff-a744-13b5e6fef561
# ╠═3e3b9f60-6791-4922-8d3a-0dbe98437060
# ╠═5379c325-69fa-4f3e-957d-644d40155997
# ╠═926a554f-e587-40fa-99b5-fad6171fd551
# ╟─6045bf4a-e6ac-4411-bd0e-171ee51f3c97
# ╠═a82986c7-007e-4b6c-9afa-530dc062c5ce
# ╠═63b2e16e-5614-472d-88d5-09ef17c6aeaa
# ╟─d2350cc9-b4a1-4ab2-ab0a-39c0dad37681
# ╠═1fbcf167-e742-48a7-a939-445ac0798550
# ╠═1b62b57e-d6c8-4fa8-8d4a-f7e62e5b795f
# ╟─5d443572-9d5e-4495-867f-e908d9b56716
# ╠═4a1aa64a-96da-4ecc-b98d-9623da487b9f
# ╠═e085b47a-a8a3-4620-9fce-f031c837e33b
# ╠═fb36dabc-4aa0-40dd-ab3b-c7b97f34ab37
# ╠═039bdb68-c4ce-4a42-ab4a-8f91c5ac29d6
# ╠═317c15e3-dc43-494a-8ad6-377b1d7a9ff8
# ╠═85a067b0-24e1-47e6-87e7-6850f981cf04
# ╠═bdf6e2d5-30df-4bb5-a755-b71f29d99b41
# ╠═bbe3111b-18b3-4832-b636-621567a19632
# ╠═980eed0d-3ba6-48a1-b31f-8a1446c96444
# ╠═e19bf203-1eac-457f-bd61-1d9bfff9aa88
# ╠═6aa9e5de-9daf-4b9a-9078-b909a9fc5e03
# ╠═18ebd13e-88cd-453c-8690-1efc55b35900
