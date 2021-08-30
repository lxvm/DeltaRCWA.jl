nbdir = ENV["HOME"] * "/.julia/dev/DeltaRCWA/notebooks/"
# begin
#     import Pkg
#     Pkg.activate(nbdir)
# end

using Revise
using DeltaRCWA
using FFTW
using LinearAlgebra
using Plots
include("$(nbdir)NystromMethodQP.jl");

struct ComplexExpSheet{T} <: RCWASheet{T, 1}
	θ::T
	θᵗ::T
	k::T
	d::T
	L::T
end

### Define how to convert between η/μ and conductivity matrix conventions
function DeltaRCWA.σₑˣˣ(sheet::ComplexExpSheet, x⃗)
	2 ./ (-sin(sheet.θ)*[1e-8+ 1-exp(1im*sheet.k*sheet.d*e[1]) for e in Iterators.product(x⃗...)])
end
function DeltaRCWA.σₘʸʸ(sheet::ComplexExpSheet, x⃗)
	-2sin(sheet.θ)*[1+exp(1im*sheet.k*sheet.d*e[1]) for e in Iterators.product(x⃗...)]
end

function compare_DeltaRCWA_BoundaryIntegral(θᵗ)
# constants
k = 10.0   # wavenumber k²=ω²μ₀ϵ₀
λ = 2*pi/k # wavelength 
θ = -π/2.0 # incidence angle (measured with respect to the x axis)
α = k*cos(θ)
β = k*sin(θ)
uInc(x,y)= exp(im*α*x+im*β*y)  # incident planewave
# θᵗ = -π/3;    # transmitted field angle
d  = cos(θᵗ)-cos(θ); 
# Defining susceptibilities
η(x) = -sin(θ).*(1+exp.(im*k*d*x));
μ(x) = 1e-8-sin(θ).*(1-exp.(im*k*d*x));
L = 2*(2*π)/(k*abs(d));  # Unit cell width

# Carlos' solver
# Defining geometry
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

# Solving the boundary integral equations
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

lims = [-L L -L L]
h=0.025
x = lims[1]:h:lims[2]
y = lims[3]:h:lims[4]

ver⁺ = [L 0;-L 0;-L 100;L 100]
Ω⁺ = isin(x,y,ver⁺)

ver⁻ = [L 0;-L 0;-L -100;L -100]
Ω⁻ =isin(x,y,ver⁻)
UInc = reshape(uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx)
@assert UInc ≈ reshape(uInc.(Ω⁻.pts[:,1],Ω⁻.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx)


pot⁺ = potentialsSkeleton(Ω⁺.pts,Γ,k,L,α,NC);
U⁺scat = reshape(pot⁺.SL*φ⁺,Ω⁺.Ny,Ω⁺.Nx)
U⁺   =  Ω⁺.In .* (U⁺scat .+ UInc)

pot⁻= potentialsSkeleton(Ω⁻.pts,Γ,k,L,α,NC);
U⁻scat = reshape(pot⁻.SL*φ⁻,Ω⁻.Ny,Ω⁻.Nx)
U⁻   =  Ω⁻.In .* (U⁻scat .+ UInc)


pltReal =  heatmap(x,y,real.(U⁺),color=:RdBu,clim=(-1.0,1.0))
pltReal =  heatmap!(x,y,real.(U⁻),color=:RdBu,clim=(-1.0,1.0))
pltReal = plot!(pltReal,legend=false,aspect_ratio=:equal,frame=:box)
display(pltReal)
pltImag =  heatmap(x,y,imag.(U⁺),color=:RdBu,clim=(-1.0,1.0))
pltImag =  heatmap!(x,y,imag.(U⁻),color=:RdBu,clim=(-1.0,1.0))
pltImag = plot!(pltImag,legend=false,aspect_ratio=:equal,frame=:box)
display(pltImag)

# heatmap(x,y,abs2.(U⁻),color=:RdBu,clim=(-1.0,1.0))
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

display(plot(sol, method=:fft, combine=false,
	# mask=Bool[n==argmax(abs.(sol.O₂)) ? false : true for n in 1:length(nvec)]
))


# wrangle Luke's solution into DeltaRCWA for plotting
Lukemodes = PlanewaveModes(k, Vacuum(), ((length(nvec), L), ), prob.modes.x⃗, (kₓ, ), βz)
Lukesol = DeltaRCWA.DeltaRCWASolution(Lukemodes, pol, u_n, u_p, u_out_n, u_out_p)
display(plot(Lukesol))


# error of scattered amplitudes in port 1
norm(u_out_n - sol.O₁)

# error of scattered amplitudes in port 2
norm(u_out_p - sol.O₂)

S ≈ smatrix(sheet, prob.modes, pol)



# find the subsets of the x, y points in each domain (Ω⁺, Ω⁻)
xisinΩ⁺ = [1.0 in Ω⁺.In[:, i] ? true : false for i in eachindex(x)]
yisinΩ⁺ = [1.0 in Ω⁺.In[i, :] ? true : false for i in eachindex(y)]
xinΩ⁺ = x[findfirst(xisinΩ⁺):findlast(xisinΩ⁺)]
yinΩ⁺ = y[findfirst(yisinΩ⁺):findlast(yisinΩ⁺)]
xisinΩ⁻ = [1.0 in Ω⁻.In[:, i] ? true : false for i in eachindex(x)]
yisinΩ⁻ = [1.0 in Ω⁻.In[i, :] ? true : false for i in eachindex(y)]
xinΩ⁻ = x[findfirst(xisinΩ⁻):findlast(xisinΩ⁻)]
yinΩ⁻ = y[findfirst(yisinΩ⁻):findlast(yisinΩ⁻)]

k⃗ = sol.modes.k⃗[1]
kz = prob.modes.kz
# Luke's equivalent
# k⃗ = kₓ
# kz = βz

mask = fill(true, size(kz))
mask = Bool[n==argmax(abs.(sol.O₂)) ? false : true for n in 1:length(nvec)]
O₂ = sol.O₂ .* mask
I₂ = sol.I₂ .* mask
I₁ = sol.I₁ .* mask
O₁ = sol.O₁ .* mask
# To compare with Luke, use these
# O₂ = Lukesol.O₂ .* mask
# I₂ = Lukesol.I₂ .* mask
# I₁ = Lukesol.I₁ .* mask
# O₁ = Lukesol.O₁ .* mask

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


# Carlos: extract the fields in the upper domain
Ω⁺mask = isfinite.(Ω⁺.In)
CI₂ = reshape(UInc[Ω⁺mask], length.((yinΩ⁺, xinΩ⁺)))
CO₂ = reshape(U⁺[Ω⁺mask] .- UInc[Ω⁺mask], length.((yinΩ⁺, xinΩ⁺)))


# Carlos: extract the fields in the lower domain
Ω⁻mask = isfinite.(Ω⁻.In)
CI₁ = reshape(UInc[Ω⁻mask] .- UInc[Ω⁻mask], length.((yinΩ⁻, xinΩ⁻)))
CO₁ = reshape(U⁻[Ω⁻mask], length.((yinΩ⁻, xinΩ⁻)))

for part in (real, imag, abs)
clim = (-1, 1)
# Incident, upper
display(heatmap(xinΩ⁺, yinΩ⁺, part.(AI₂ .- CI₂), clim=clim, xguide="x", yguide="z", title="$(part)(error incident field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
display(histogram(reshape(part.(AI₂ .- CI₂), :), normalize=true, title="density $(part)(error incident field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
# scattered, upper
display(heatmap(xinΩ⁺, yinΩ⁺, part.(AO₂ .- CO₂), clim=clim, xguide="x", yguide="z", title="$(part)(error scattered field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
display(histogram(reshape(part.(AO₂ .- CO₂), :), normalize=true, title="density $(part)(error scattered field above sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
# Incident, lower
display(heatmap(xinΩ⁻, yinΩ⁻, part.(AI₁ .- CI₁), clim=clim, xguide="x", yguide="z", title="$(part)(error incident field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
display(histogram(reshape(part.(AI₁ .- CI₁), :), normalize=true, title="density $(part)(error incident field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
# scattered, lower
display(heatmap(xinΩ⁻, yinΩ⁻, part.(AO₁ .- CO₁), clim=clim, xguide="x", yguide="z", title="$(part)(error scattered field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
display(histogram(reshape(part.(AO₁ .- CO₁), :), normalize=true, title="density $(part)(error scattered field under sheet), θᵗ=$(round(180 *θᵗ/π; digits=2))ᵒ"))
end
end