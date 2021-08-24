### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 5c283082-518b-45c6-907a-e80945b1a4e4
begin
	using Pkg
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

For references on the test choice see
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

# The parameter M defines the number of boundary points to used in the discretization of the BIE
M  = Int(round(max(20,20/λ)))

# Γ contains all the information of the discretized curve that is needed for the numerical solution of the problem
Γ  = SkeletonCctor(ver,info,[M])
end;

# ╔═╡ 1710cfd1-89c3-4bf8-b825-e76e7c898d74
begin
xPlot = -L/2:0.001:L/2;

paramPlot = plot(xPlot,real.(η.(xPlot)))

paramPlot = plot!(xPlot,imag.(η.(xPlot))) 

paramPlot = plot!(xPlot,real.(μ.(xPlot))) 

paramPlot = plot!(xPlot,imag.(μ.(xPlot))) 

plot(paramPlot,lw=3,linestyle=[:solid :dash :solid :dash],label = ["Re η" "Im η" "Re μ" "Im μ"],frame=:box)
xlabel!("x");ylabel!("η, μ")
end

# ╔═╡ 17f4f76d-fadd-4f0b-84c0-8d5dd8582176
f1(x) = -2*im*k*η(x)*exp(im*α*x)

# ╔═╡ 7d2b7703-8001-4308-9ca8-bf736d40099b
f2(x) = -2*im*β*exp(im*α*x)

# ╔═╡ b08df8b1-de6f-4b71-a650-48e483d91a3e
# Solving the boundary integral equations
begin
method = "MK" # Selection of the Nystrom method to be used

NC  = 5;      # Number of unit cells considered in the windowed sum (A = NC*L)

mat = matricesSkeleton(1,Γ,k,L,α,NC,method) # Construction of the Nystrom matrices
    
x1 = Γ.parts[1].x[:,1]; # x-coordinate of the discretization points on Γ
    
Id = Matrix(I,Γ.Npts,Γ.Npts); # Indentity matrix

A⁺ = -0.5*Id+im*k*Diagonal(η.(x1))*mat.S  # System matrix for Ω⁺
b⁺ = f1.(x1) # boundary data for Ω⁺
ψ⁺ = A⁺\b⁺; 

A⁻ = -0.5*Id+im*k*Diagonal(μ.(x1))*mat.S # System matrix for Ω⁻
b⁻ =  f2.(x1) # boundary data for Ω⁻
ψ⁻ = A⁻\b⁻; 

φ⁺ = 0.5*(ψ⁺+ψ⁻);
φ⁻ = 0.5*(ψ⁺-ψ⁻);
end 

# ╔═╡ 1aba587b-475a-49cb-aebe-06a7e1929c5a
lims = [-L L -L L]

# ╔═╡ ed7297fb-44ba-4efe-b72c-1abb80bb7484
h=0.025

# ╔═╡ 0fd9999a-e379-46a0-8403-0fbae721b19b
x = lims[1]:h:lims[2]

# ╔═╡ 6fdbb254-c11b-4688-be59-6e909790d620
y = lims[3]:h:lims[4]

# ╔═╡ ec2d9dea-e4d7-4703-8422-f4d9c4aec3fd
begin
ver⁺ = [L 0;-L 0;-L 100;L 100]; 
Ω⁺ = isin(x,y,ver⁺);

ver⁻ = [L 0;-L 0;-L -100;L -100]; 
Ω⁻ =isin(x,y,ver⁻);
end;

# ╔═╡ 1fd8b4a3-53ef-4752-b817-bda9978c8a8e
plot(real.(ifft(uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]))))

# ╔═╡ 2af29afd-37ea-4abf-b0c4-737cfd96e101
begin
pot⁺ = potentialsSkeleton(Ω⁺.pts,Γ,k,L,α,NC) 
U⁺   =  Ω⁺.In.*reshape(pot⁺.SL*φ⁺ + uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx);

pot⁻= potentialsSkeleton(Ω⁻.pts,Γ,k,L,α,NC) 
U⁻  =  Ω⁻.In.*reshape(pot⁻.SL*φ⁻ + uInc.(Ω⁻.pts[:,1],Ω⁻.pts[:,2]),Ω⁻.Ny,Ω⁻.Nx);
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
## Comparison with DeltaRCWA
"

# ╔═╡ 6045bf4a-e6ac-4411-bd0e-171ee51f3c97
md"
## Luke's program

Extra test: add a trivial delta layer that should transmit everything and check that 
it doesn't add a phase normalization
"

# ╔═╡ a82986c7-007e-4b6c-9afa-530dc062c5ce
begin
# Some parameters repeated from Carlos
# k = 10.0;
# λ = 2*pi/k;
# θ = -π/2.0; # incidence angle (measured with respect to the x axis)
# α = k*cos(θ);
# β = k*sin(θ);
uinc(x,y)= @. exp(1im*α*x+1im*β*y);  # incident planewave

# θᵗ = -π/8;    # transmitted field angle
# d  = cos(θᵗ)-cos(θ); 
# L = 2*(2*π)/(k*abs(d));  # Unit cell width
# M₀(x) =  [1e-8 for i in x] #-sin(θ)*(1-exp(1im*k*d*x));
M₀(x) = @. -sin(θ)*(1+exp(1im*k*d*x));
N₀(x) = @. -sin(θ)*(1-exp(1im*k*d*x));
# N₀(x) =  [-10000000 for i in x] #-sin(θ)*(1-exp(1im*k*d*x));

nvec = 0:99;
dx = L/length(nvec)
xvec = [n*dx-L/2 for n in nvec]
kₓ = 2*pi*fftfreq(length(nvec), 1/dx)
βz = @. sqrt(Complex(k^2 - kₓ^2))

u_p = fft(uinc(xvec, 0))/length(nvec)
u_n = zeros(length(nvec))
MM = Matrix(Diagonal(M₀(xvec)))
# MM = Matrix(Diagonal(im * ones(size(xvec, 1)))
	
NN = Matrix(Diagonal(N₀(xvec)))
# display(N)


A = [-Diagonal(βz)-k*ifft(fft(MM, 2), 1)    -Diagonal(βz)-k*ifft(fft(MM, 2), 1);
    -Diagonal(βz)-k*ifft(fft(NN, 2), 1)     Diagonal(βz)+k*ifft(fft(NN, 2), 1)]
B = [-Diagonal(βz)+k*ifft(fft(MM, 2), 1)    -Diagonal(βz)+k*ifft(fft(MM, 2), 1);
    -Diagonal(βz)+k*ifft(fft(NN, 2), 1)     Diagonal(βz)-k*ifft(fft(NN, 2), 1)]
S = A\B
# display(S)

u_in = [u_p; u_n]
# smat_star(smat_star(smatrix(Air(L), kxy, ω₀), S), smatrix(Air(L), kxy, ω₀))
u_out = S*u_in
u_out_p = u_out[1:length(nvec)]
u_out_n = u_out[length(nvec)+1:2*length(nvec)]
end

# ╔═╡ 07b0e662-b295-46ff-a744-13b5e6fef561
begin
### DeltaRCWA solver
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
end;

# ╔═╡ 382b1cb8-45ab-4098-b7c0-f7826eba592d
prob = DeltaRCWAProblem(sheet, dims, ω, pol, zeros(length(xvec)), ifft(uInc.(xvec, 0)))

# ╔═╡ 5379c325-69fa-4f3e-957d-644d40155997
any(iszero.(prob.modes.kz))

# ╔═╡ cd9bf171-5809-4f89-8fd6-0216e9bcec01
sol = solve(prob)

# ╔═╡ 3e3b9f60-6791-4922-8d3a-0dbe98437060
plot(sol, method=:fft, combine=false)

# ╔═╡ 7246eb23-1c2a-4bde-abd5-1b996977bb5d
abs2.(sol.I₂)

# ╔═╡ bf690896-e606-4e6e-a950-1ebc96a23ddb
# plot(abs2.(sol.O₂))
plot(real.(log.(sol.I₂)))

# ╔═╡ 37cba2fa-e61e-4172-8840-28978bf9d16b
argmax(imag.(prob.modes.kz))

# ╔═╡ 5eeb32ad-07e7-43c8-b78a-d23177eb5fc4
begin
LukexPlot = -L/2:0.001:L/2;

LukeparamPlot = plot(LukexPlot,real.(M₀(LukexPlot)))

LukeparamPlot = plot!(LukexPlot,imag.(M₀(LukexPlot))) 

LukeparamPlot = plot!(LukexPlot,real.(N₀(LukexPlot))) 

LukeparamPlot = plot!(LukexPlot,imag.(N₀(LukexPlot))) 

plot(LukeparamPlot,lw=3,linestyle=[:solid :dash :solid :dash],label = ["Re M" "Im M" "Re N" "Im N"],frame=:box)
xlabel!("x");ylabel!("M, N")
end

# ╔═╡ ebef37ce-d77e-4fc0-84d2-d147c6f9a0b6
begin
E1_in(x,z) = sum([u_p[n+1]*exp(-1im*kₓ[n+1]*x)*exp(1im*βz[n+1]*z) for n in nvec])
E1_out(x,z) = sum([u_out_p[n+1]*exp(-1im*kₓ[n+1]*x)*exp(1im*βz[n+1]*z) for n in nvec])
E1(x, z) = E1_in(x, z) + E1_out(x, z)
E2(x,z) = sum([u_out_n[n+1]*exp(-1im*kₓ[n+1]*x)*exp(-1im*βz[n+1]*z) for n in nvec])

xview = range(-L, L, length=size(U⁺, 1)) #-L:0.02:L 
zview = range(0, L, length=size(U⁻, 2) ÷ 2) #0:0.01:L
zview2 = range(-L, 0, length=size(U⁺, 2) ÷ 2) #-L:0.01:0
Emat_in = zeros(Complex{Float64}, 0)
Emat_out = zeros(Complex{Float64}, 0)
Emat = zeros(Complex{Float64}, 0)
Emat2 = zeros(Complex{Float64}, 0)
for (z1, z2) in zip(zview, zview2)
    append!(Emat_in, E1_in.(xview, z1))
    append!(Emat_out, E1_out.(xview, z1))
    append!(Emat2, E2.(xview, z2))
end
Emat = Emat_in + Emat_out
Emat = reshape(Emat,length(xview),:)
Emat2 = reshape(Emat2,length(xview),:)
end;

# ╔═╡ 68b8af24-098a-4a18-be7a-2245a787d325
begin
plt1 = heatmap(xview, zview2, color=:RdBu, transpose(real.(Emat2)))
plt1 = heatmap!(xview, zview, color=:RdBu, transpose(real.(Emat)))
plot(plt1)
end

# ╔═╡ 161f163b-7b4b-443f-98ce-b8436f045d5b
begin
plt2 = heatmap(xview, zview2, color=:RdBu, transpose(imag.(Emat2)))
plt2 = heatmap!(xview, zview, color=:RdBu, transpose(imag.(Emat)))
plot(plt2)
end

# ╔═╡ a9b675ab-8346-4ed4-acb1-6fd294c887a4
begin
plt3 = heatmap(xview, zview2, color=:RdBu, transpose(abs2.(Emat2)))
plt3 = heatmap!(xview, zview, color=:RdBu, transpose(abs2.(Emat)))
plot(plt3)
end

# ╔═╡ 7616b15c-e4ba-4276-aaf1-6db59b3e971e
md"
## Comparisons

First we should extract the difference in normalization between the two models.
This can be done by comparing the amplitudes of the incident waves, since propagation
in free space is linear.
We can just take the maximum of the incident wave
"

# ╔═╡ 9d054936-196d-4019-87ea-4e1005776b4d
Lukemax = maximum(real.(Emat_in))

# ╔═╡ ca0e851b-c86e-4bad-9475-2d5e2724a3dc
Carlosmax = maximum(real.(uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2])))

# ╔═╡ 843a7cf5-293e-471e-bfdc-a7b2c3e6a978
md"
Check for compatibility of element-wise field comparison.
May need to interpolate between the grids
"

# ╔═╡ 7815c22d-2351-48f4-9ba5-57654bc802c3
[size(e) for e in [U⁺, U⁻, Emat, Emat2, xview, zview, zview2]]

# ╔═╡ 6b782b97-ef23-4831-9871-5b4c46d03e91
begin
fielderrorplot = heatmap(xview, zview2, color=:RdBu, abs.(transpose(Emat2)))
# fielderrorplot = heatmap!(xview, zview, color=:RdBu, abs.(U⁺ - transpose(Emat_in)))
plot(fielderrorplot)
end
# find a global phase difference

# ╔═╡ 0c3fd6a9-d494-4d7c-bdb7-29823b83c798
begin
	phaseplot = heatmap(xview, zview, angle.(transpose(Emat_in)))
	plot(phaseplot)
end

# ╔═╡ 0d471a35-7b73-4520-b899-324343b371cd
begin
	pphaseplot = heatmap(x, y, angle.(reshape(uInc.(Ω⁻.pts[:,1],Ω⁻.pts[:,2]),Ω⁻.Ny,Ω⁻.Nx)))
	plot(pphaseplot)
end

# ╔═╡ b27d8534-0173-4831-8f72-1c9b93eda72e
begin
	ppphaseplot = heatmap(x, y, angle.(reshape(uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx)))
	plot(ppphaseplot)
end


# ╔═╡ 15d11b59-0a4d-4f18-8704-3277142f8b59
size(reshape(uInc.(Ω⁻.pts[:,1],Ω⁻.pts[:,2]),Ω⁻.Ny,Ω⁻.Nx))
# Ω⁻.pts[:,1]
# Ω⁻.pts[:,2]

# ╔═╡ 3199a72b-6c1d-4479-96f2-44ec877f8a4f
# show the difference in phase
size(reshape(uInc.(Ω⁻.pts[:,1],Ω⁻.pts[:,2]),Ω⁻.Ny,Ω⁻.Nx)[1, :])

# ╔═╡ 65008c5a-f0d1-471e-81d2-037e0c577ecf
size(transpose(Emat_in))

# ╔═╡ Cell order:
# ╠═a8253062-df3c-11eb-26e1-0dcff18bf886
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
# ╠═17f4f76d-fadd-4f0b-84c0-8d5dd8582176
# ╠═7d2b7703-8001-4308-9ca8-bf736d40099b
# ╠═1aba587b-475a-49cb-aebe-06a7e1929c5a
# ╠═ed7297fb-44ba-4efe-b72c-1abb80bb7484
# ╠═0fd9999a-e379-46a0-8403-0fbae721b19b
# ╠═6fdbb254-c11b-4688-be59-6e909790d620
# ╠═ec2d9dea-e4d7-4703-8422-f4d9c4aec3fd
# ╠═1fd8b4a3-53ef-4752-b817-bda9978c8a8e
# ╠═2af29afd-37ea-4abf-b0c4-737cfd96e101
# ╠═a39dd1aa-868d-4819-877f-ccd4b856cb4c
# ╠═a022ef8b-03cc-42cc-ac9e-4da96e3edecb
# ╠═5379c325-69fa-4f3e-957d-644d40155997
# ╠═66a8554e-1b71-4776-a903-c9ee54f1d64b
# ╟─78e6c13a-ba9a-4b5f-8ae6-51ebeb7f098e
# ╠═07b0e662-b295-46ff-a744-13b5e6fef561
# ╠═382b1cb8-45ab-4098-b7c0-f7826eba592d
# ╠═cd9bf171-5809-4f89-8fd6-0216e9bcec01
# ╠═3e3b9f60-6791-4922-8d3a-0dbe98437060
# ╠═7246eb23-1c2a-4bde-abd5-1b996977bb5d
# ╠═37cba2fa-e61e-4172-8840-28978bf9d16b
# ╠═bf690896-e606-4e6e-a950-1ebc96a23ddb
# ╟─6045bf4a-e6ac-4411-bd0e-171ee51f3c97
# ╠═a82986c7-007e-4b6c-9afa-530dc062c5ce
# ╠═5eeb32ad-07e7-43c8-b78a-d23177eb5fc4
# ╠═ebef37ce-d77e-4fc0-84d2-d147c6f9a0b6
# ╠═68b8af24-098a-4a18-be7a-2245a787d325
# ╠═161f163b-7b4b-443f-98ce-b8436f045d5b
# ╠═a9b675ab-8346-4ed4-acb1-6fd294c887a4
# ╠═7616b15c-e4ba-4276-aaf1-6db59b3e971e
# ╠═9d054936-196d-4019-87ea-4e1005776b4d
# ╠═ca0e851b-c86e-4bad-9475-2d5e2724a3dc
# ╠═843a7cf5-293e-471e-bfdc-a7b2c3e6a978
# ╠═7815c22d-2351-48f4-9ba5-57654bc802c3
# ╠═6b782b97-ef23-4831-9871-5b4c46d03e91
# ╠═0c3fd6a9-d494-4d7c-bdb7-29823b83c798
# ╠═0d471a35-7b73-4520-b899-324343b371cd
# ╠═b27d8534-0173-4831-8f72-1c9b93eda72e
# ╠═15d11b59-0a4d-4f18-8704-3277142f8b59
# ╠═3199a72b-6c1d-4479-96f2-44ec877f8a4f
# ╠═65008c5a-f0d1-471e-81d2-037e0c577ecf
