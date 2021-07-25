### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 5c283082-518b-45c6-907a-e80945b1a4e4
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 898fa516-591d-4831-8c65-f0e758922d15
using DeltaRCWA

# ╔═╡ 8f32367b-b19b-4c18-80ff-5edaabcfb0d9
include("NystromMethodQP.jl")

# ╔═╡ a8253062-df3c-11eb-26e1-0dcff18bf886
md"
# Testing multilayer films
"

# ╔═╡ e944920d-7aa0-4fdc-ad47-4c49821eceeb
# unit of length
L₀ = 1

# ╔═╡ ca120085-b234-48aa-9693-47171a799ea2
md"
## Setup layers

I will assemble a structure of two air layers and propagate waves through them.
The structure is 2D and periodic.
"

# ╔═╡ 7cc184e4-9f8d-4f2e-82cf-f5872b115102
# Depth of layer 1
d₁ = 1 * L₀

# ╔═╡ 618bd4c7-a089-4d25-9491-16433c5479ae
# First air layer
# air₁ = Air(d₁)
air₁ = UniformLayer(d₁, 1, 1)

# ╔═╡ fc3390a9-caef-4602-b9fe-fe83cbe17cad
# Depth of layer 2
d₂ = 1 * L₀

# ╔═╡ f0580da9-52db-4c40-8706-b537870ab1dc
air₂ = UniformLayer(1, 1, 1)

# ╔═╡ 4eb88559-8bbc-49f8-8d46-a4485892f57f
md"
## Define spatial x, y lattices

The normally incident wave will align with the z axis, and the x axis will be
periodic, whereas the y axis will be invariant.

However, I will still write the code for a 2D photonic crystal, with a trivial y axis
"

# ╔═╡ 159bdc89-6e78-435d-81e3-3cb2b77312e4
Nxmodes = 10

# ╔═╡ 153b9330-b0c4-4b84-9a26-aff06089b334
xmodenumbers = 0:(Nxmodes - 1)

# ╔═╡ 02c7d3b2-f941-4973-8339-e5eb8bda0f19
# x lattice vector (period of unit cell)
ax = 1 * L₀

# ╔═╡ d8e567cb-90b8-4a36-902a-a980c684f728
# spacing between the x grid points 
δx = ax / Nxmodes

# ╔═╡ c2503f86-3083-4ebf-8a17-350bb1748e97
# x grid points
begin
	xgrid = cumsum([δx for i in 0:(Nxmodes - 2)])
	pushfirst!(xgrid, 0)
end

# ╔═╡ 25396236-c4ff-4161-8d1b-edb07174de2c
# 1 mode means translation invariant
Nymodes = 1

# ╔═╡ dee4e053-cc7e-41ad-ba24-2a641131ce38
ymodenumbers = 0:(Nymodes - 1)

# ╔═╡ 9cf2ab2a-e9ab-4194-bdf3-5db2c4939a03
# y lattice vector (period of unit cell)
ay = 1 * L₀

# ╔═╡ 4ff1f5e6-1cce-4da1-b2ef-99e491ea5719
# spacing between the y grid points 
δy = ay / (Nymodes - 1) 

# ╔═╡ ade9e849-35ed-4943-9fd8-49f32d5fea61
# y grid points
begin
	ygrid = cumsum([δy for i in 0:(Nymodes - 2)])
	pushfirst!(ygrid, 0)
end

# ╔═╡ e89f082d-8cdb-4bc9-8d88-03e57c1fc11a
md"
## Define reciprocal x, y lattices

Setup the normal modes in the x direction, due to periodicity

From 'Photonic Crystals', we are doing the decomposition

$\boldsymbol H_{(n, k_z, \boldsymbol k_\parallel)} (\boldsymbol r) = 
e^{i \boldsymbol k_\parallel \cdot \boldsymbol \rho} e^{i k_z z}
\boldsymbol u_{(n, k_z, \boldsymbol k_\parallel)} (\boldsymbol \rho)$
where $\boldsymbol k_\parallel = (k_x, k_y)$ and $\boldsymbol \rho = (x, y)$.
Now this is totally valid for 1D and 2D photonic crystals, where the z direction is 
translation invariant.
For our problem, the z direction will be all air except for sets of Lebesgue measure
zero for which scattering planes exist.
We need to match boundary conditions at these scattering planes, which may require
several modes $k_z$ in addition to the modes $\boldsymbol k_\parallel$ along the
surface.

We are using only a 1D photonic crystal for this example, so we neglect $k_y$.
"

# ╔═╡ 91b692b2-bfff-4119-84ec-c3f6dc2e9120
# Set the modes for the x periodic expansion
kx = fftfreq(Nxmodes, Nxmodes / ax)

# ╔═╡ 6e21e9a8-c7f0-474c-9730-2a0d4c768fb9
# Set the modes for the y periodic expansion
ky = fftfreq(Nymodes, Nymodes / ay)

# ╔═╡ 63b477fd-38c5-4532-a752-916a5f08bee7
# Define the wavevectors for the periodic modes parallel to scattering planes
begin
	Nxymodes = Nxmodes * Nymodes
	kxy = zeros(2, Nxymodes)
	local i = 0
	for ey in ky
		for ex in kx
			i += 1
			kxy[1, i] = ex
			kxy[2, i] = ey
		end
	end
	# so kxy[:, i] gives the ith wavevector
end

# ╔═╡ 53b9bb4a-96b9-40f5-85f5-42d6e0302703
md"
## Setup the incident modes on the system

In addition we choose a total frequency for these modes, ω₀, with wavelength
λ₀.
"

# ╔═╡ 815d975c-0ad4-491d-bd2e-8f0053ffa29b
λ₀ = 1.2 * L₀

# ╔═╡ fbd43517-8d23-4d53-8919-584a4b45817c
ω₀ = 2 * pi / λ₀

# ╔═╡ 4e78b99e-8084-4ed9-886d-f9bfd3190640
# x mode amplitudes incident on region 1
# set all even modes to 1
cxin₁ = [mod(n, 2) == 0 ? 1 : 0 for n in xmodenumbers]

# ╔═╡ 801cf9ef-fe46-480a-9124-601fa03f62da
# x mode amplitudes incident on region 2
# set all even modes to 1
cxin₂ = [mod(n, 2) == 0 ? 1 : 0  for n in xmodenumbers]

# ╔═╡ fd1f4984-a0d0-498d-832f-1f884acf0265
# y mode amplitudes incident on region 1
# set all modes to 1
cyin₁ = [1 for n in ymodenumbers]

# ╔═╡ aaa1d158-59b3-4c72-9306-fb7bcef3a09e
# y mode amplitudes incident on region 1
# set all modes to 1
cyin₂ = [1 for n in ymodenumbers]

# ╔═╡ 65ffffaa-6511-49df-b0ed-0353c52113d6
# get total incident mode matrix on region 1
cxyin₁ = kron(cyin₁, cxin₁)

# ╔═╡ 3108c320-3f17-4e0b-858b-b9e25c11f02e
# get total incident mode matrix on region 2
cxyin₂ = kron(cyin₂, cxin₂)

# ╔═╡ 484f4b6a-8ccb-4519-9d04-fb0b505fc13a
cxyin = mortar([cxyin₁, cxyin₂])

# ╔═╡ f47a4fcb-fa16-48ab-850c-2a3de6bfeac7
md"
## Propagate waves
"

# ╔═╡ 6d4dd4c0-8e0c-4825-aa64-3a82ced17cd3
sair₁ = smatrix(air₁, (kx, ky), ω₀) 
# == smatrix(air₁, reshape([i for i in kx], 1, size(kx, 1)), ω₀)

# ╔═╡ 90d06519-9526-48c0-8560-116b7c112997
reshape(([e[1] + e[2] for e in Iterators.product((kx, ky)...)]), prod(length, (kx, ky)))

# ╔═╡ 9ba33324-8561-4413-8343-f9b52601d788
size(kx)

# ╔═╡ 0342e0f7-e654-41c2-9a56-b5e99d59ca27
sair₂ = smatrix(air₂, (kx, ky), ω₀) 

# ╔═╡ 36a24a19-bbdc-466e-8905-da41493c5f8d
sair₁₂ = smat_star(sair₁, sair₂)

# ╔═╡ 7fdb6893-010e-4a6a-9636-d59e966e2fac
sair₁₂[Block(1, 2)]

# ╔═╡ 3db65194-8647-41ff-8d30-8c79638b3753
# Alternative construction
airstack = ScatteringStack([air₁, air₂])

# ╔═╡ 06a9e35c-0400-4875-9ea7-86927ff288fa
smatrix(airstack, (kx, ky), ω₀) == sair₁₂ 

# ╔═╡ 315bffb4-ef0f-4dcf-8b0c-31fbd2f0a795
(sair₁₂ * cxyin)[Block(1)]

# ╔═╡ 7229e888-dc22-47d2-8dd5-db06bd5ae476
md"
### Check unitarity of scattering matrices

I should really be checking unitarity with respect to the star product, but this 
example is with diagonal scattering matrices and you can't tell the difference
"

# ╔═╡ 9a7c15d1-6ac9-4e6b-851f-1ef502507254
[norm(M'M - I) / norm(M'M) for M in [sair₁, sair₂, sair₁₂]]

# ╔═╡ 22081b08-f086-4a64-b349-b6938c1c4da1
md"
## Build scattering matrix for delta layer

The following explanation follow Luke's notes

We need functions that get the M and N matrices for each polarization
From Luke's notes:
- $k = \sqrt{\epsilon_0 \mu_0}\omega, \eta = \sqrt{\mu_0 / \epsilon_0}$
- For TE polarization $M = \frac{\eta}{2Z}, N = 2 Y \eta$
- For TM polarization $M = \frac{1}{2Y\eta}, N = \frac{2 Z}{\eta}$
Question! if $Y=Z^{-1}$ then is this a matrix inversion or element-wise reciprocation
Answer: Matrix inversion (but which is element-wise reciprocation for diagonal matrices)

In each case, if the field component of interest (The one which spans only one basis
vector in each TE and TM polarization) is labelled by $\boldsymbol c$ then the
sheet transition condition requires

$[[\partial_z \boldsymbol c]] = -i k M \{\!\{\boldsymbol c\}\!\}$
$\{\!\{\partial_z \boldsymbol c\}\!\} = -i k N [[\boldsymbol c]]$

where $\{\!\{ \boldsymbol c \}\!\} = \boldsymbol c_1 + \boldsymbol c_2$
and  $[[ \boldsymbol c ]] = \boldsymbol c_1 - \boldsymbol c_2$.
Here the subscript denotes the port of the scatterer at which the field is evalutated.

We may decompose the field into left and right-propagating components
(alternatively incident and scattered components) by
$\boldsymbol c_i = \boldsymbol c_i^+ + \boldsymbol c_i^-$.
Considering these conventions, we can rewrite the boundary conditions as


$\partial_z
(\boldsymbol c_1^+ + \boldsymbol c_1^- - \boldsymbol c_2^+ - \boldsymbol c_2^-)
= -i k M
(\boldsymbol c_1^+ + \boldsymbol c_1^- + \boldsymbol c_2^+ + \boldsymbol c_2^-)$
$\partial_z 
(\boldsymbol c_1^+ + \boldsymbol c_1^- + \boldsymbol c_2^+ + \boldsymbol c_2^-)
= -i k N
(\boldsymbol c_1^+ + \boldsymbol c_1^- - \boldsymbol c_2^+ - \boldsymbol c_2^-)$

which can be rearranged into this linear system

$\begin{pmatrix}
\partial_z + ikM  & -\partial_z + ikM
\\
\partial_z + ikN & \partial_z - ikN
\end{pmatrix}
\begin{pmatrix}
\boldsymbol c_1^-
\\
\boldsymbol c_2^+
\end{pmatrix}
= 
-\begin{pmatrix}
\partial_z + ikM  & -\partial_z + ikM
\\
\partial_z + ikN & \partial_z - ikN
\end{pmatrix}
\begin{pmatrix}
\boldsymbol c_1^+
\\
\boldsymbol c_2^-
\end{pmatrix}$

We now want to apply the plane-wave ansatz for the propagation of the fields in the
z-direction which is that $\boldsymbol c_j^\pm \propto e^{\pm i \beta_j z}$.
Hence $\partial_z \boldsymbol c_j^\pm = \pm i \beta_j \boldsymbol c_j^\pm$.
We will apply an additional transformation, that of changing the linear operators
acting on the x axis into the Fourier basis.
In total, the linear system can be rewritten as:
$\begin{pmatrix}
-i\beta_1 + ik\tilde{M}  & -i\beta_2 + ik\tilde{M}
\\
-i\beta_1 + ik\tilde{N} & i\beta_2 - ik\tilde{N}
\end{pmatrix}
\begin{pmatrix}
\tilde{\boldsymbol c}_1^-
\\
\tilde{\boldsymbol c}_2^+
\end{pmatrix}
= 
-\begin{pmatrix}
i\beta_1 + ik\tilde{M} & i\beta_2 + ik\tilde{M}
\\
i\beta_1 + ik\tilde{N} & -i\beta_2 - ik\tilde{N}
\end{pmatrix}
\begin{pmatrix}
\tilde{\boldsymbol c}_1^+
\\
\tilde{\boldsymbol c}_2^-
\end{pmatrix}$

The tensors with a tilde are in the Fourier basis and are obtained from
`M̃ M \wideutilde = ifft(fft(M, 2), 1)` for matrices.
Note that the second axis are the rows and the first are the columns.
(Hence array indexing is little-endian, since the column is rightmost)

Question! why is this it?
Answer: This is applying the Fourier operator and its adjoint to the raised and
lowered indices of the linear transformation M to get the transformed indices.
So $\tilde{M} = \mathcal F(M) = F M F^*$, where $F$ is the DFT matrix.
We can do this fast when doing a matrix-vector multiplication because the $F^*$
can be applied to the vector ($\mathcal O(n \log(n))$), then multiply by the diagonal matrix $M$ ($\mathcal O(n)$) followed by a second DFT ($\mathcal O(n \log(n))$),
which has better complexity than the $\mathcal O(n^2)$ operations we have now.

The resulting scattering matrix is obtained from inverting the left matrix
and moving it to the right side, thus expressing the inputs in terms of the outputs.

Note that the $\beta_j$ are diagonal matrices with one entry per mode.
"

# ╔═╡ 7a8ddcf3-40c4-4021-a197-8068d8cf22a2
function get_XY(layer::DeltaScatterer, xgrid, ygrid) 
	# In practice, dispatch on the concrete type of `layer`
	# to make use of the parametrization of its geometry
	error("Not implemented")
end

# ╔═╡ cc8c8ab2-d5b8-485a-a1e1-ce7e735bdf9c
function get_MN(layer::DeltaScatterer, xgrid, ygrid) 
	# In practice, dispatch on the concrete type of `layer`
	# to make use of the parametrization of its geometry
	error("Not implemented")
end

# ╔═╡ efe96bea-feee-413e-a03e-1254bc27e8d8
# use Luke's function get_all_modes for the scattering matrix

# ╔═╡ ed5050f9-8f56-475c-9abd-d661a49aceb9
md"
## Plot fields

We need to take the inverse FT of the output modes and propagate them through the
medium at the desired points (see Luke's code)
"

# ╔═╡ 1ef4e805-3dd9-4f39-a697-ddc0b4f04ad8
md"
# Testing with Carlos' solver

## Carlos' code
"

# ╔═╡ ad90ebc0-99f2-4636-bce8-12abf791a69a
# Defining constants
begin
k = 10.0;   # wavenumber k²=ω²μ₀ϵ₀

λ = 2*pi/k; # wavelength 

θ = -π/2.0; # incidence angle (measured with respect to the x axis)

α = k*cos(θ);

β = k*sin(θ);

uInc(x,y)= exp(im*α*x+im*β*y);  # incident planewave
end

# ╔═╡ fea6690d-7602-462b-85dd-4fd504918933
# Defining susceptibilities
begin
θᵗ = -π/3;    # transmitted field angle

d  = cos(θᵗ)-cos(θ); 

# η(x) = -sin(θ).*(1+exp.(im*k*d*x));
η(x) = 0#-sin(θ).*(1+exp.(im*k*d*x));

# μ(x) = -sin(θ).*(1-exp.(im*k*d*x));
μ(x) = 1000000000 #-sin(θ).*(1-exp.(im*k*d*x));

L = 2*(2*π)/(k*abs(d));  # Unit cell width
end

# ╔═╡ 1710cfd1-89c3-4bf8-b825-e76e7c898d74
begin
using Plots

xPlot = -L/2:0.001:L/2;

paramPlot = plot(xPlot,real.(η.(xPlot)))

paramPlot = plot!(xPlot,imag.(η.(xPlot))) 

paramPlot = plot!(xPlot,real.(μ.(xPlot))) 

paramPlot = plot!(xPlot,imag.(μ.(xPlot))) 

plot(paramPlot,lw=3,linestyle=[:solid :dash :solid :dash],label = ["Re η" "Im η" "Re μ" "Im μ"],frame=:box)
xlabel!("x");ylabel!("η, μ")
end

# ╔═╡ da3f4379-7b2d-4efc-a9b2-5dcccd4a7d40
# Defining geometry
begin
# Vertices that define the curve Γ 
ver = [  L/2   0;   #1 
        -L/2   0];  #2 
# In this case, Γ is just a line segment from (L/2,0) to (-L/2,0). 

# This matrix defines the orientation of the curve
        #starting node  #ending node  #domain on the right #domian on the left  
info = [1               2             1                    0];

# The parameter M defines the number of boundary points to used in the discretization of the BIE
M  = Int(round(max(20,20/λ))); 

# Γ contains all the information of the discretized curve that is needed for the numerical solution of the problem
Γ  = SkeletonCctor(ver,info,[M]);
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
end

# ╔═╡ 2af29afd-37ea-4abf-b0c4-737cfd96e101
begin
pot⁺ = potentialsSkeleton(Ω⁺.pts,Γ,k,L,α,NC) 
U⁺   =  Ω⁺.In.*reshape(pot⁺.SL*φ⁺ + uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx);

pot⁻= potentialsSkeleton(Ω⁻.pts,Γ,k,L,α,NC) 
U⁻  =  Ω⁻.In.*reshape(pot⁻.SL*φ⁻ + uInc.(Ω⁻.pts[:,1],Ω⁻.pts[:,2]),Ω⁻.Ny,Ω⁻.Nx);
end

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

# ╔═╡ 6045bf4a-e6ac-4411-bd0e-171ee51f3c97
md"
## Luke's reproduction

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
M₀(x) =  [1e-8 for i in x] #-sin(θ)*(1-exp(1im*k*d*x));
# M₀(x) = @. -sin(θ)*(1+exp(1im*k*d*x));
# N₀(x) = @. -sin(θ)*(1-exp(1im*k*d*x));
N₀(x) =  [-10000000 for i in x] #-sin(θ)*(1-exp(1im*k*d*x));

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

# ╔═╡ 95cb80ae-4979-47e1-a798-11c8f33259f5
reshape(kₓ, (size(kₓ, 1), 1))

# ╔═╡ 5eeb32ad-07e7-43c8-b78a-d23177eb5fc4
begin
LukexPlot = -L/2:0.001:L/2;

LukeparamPlot = plot(LukexPlot,real.(M₀.(LukexPlot)))

LukeparamPlot = plot!(LukexPlot,imag.(M₀.(LukexPlot))) 

LukeparamPlot = plot!(LukexPlot,real.(N₀.(LukexPlot))) 

LukeparamPlot = plot!(LukexPlot,imag.(N₀.(LukexPlot))) 

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
end

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

# ╔═╡ 992eb33a-3fcc-4d1b-bf48-3f33b8c2a1e4
md"
## Test for a trivial layer

We want to make sure that the scattering matrix for a transparent layer of delta
function thickness merely returns the identity.
"

# ╔═╡ 543c8328-ab82-435a-aca7-1f0ab573294d
struct TrivialDeltaLayer <: DeltaScatterer end

# ╔═╡ 8da39a74-6081-4bb1-8d6a-1a0d9385e189
function mysmatrix(layer::TrivialDeltaLayer, kxy::AbstractMatrix, ω₀::Real)
	Nk = size(kxy, 2)
	return mortar(reshape([
		zeros(Bool, Nk, Nk),
		Matrix(I, Nk, Nk),
		Matrix(I, Nk, Nk),
		zeros(Bool, Nk, Nk),
	], 2, 2))		
end

# ╔═╡ 521f8552-1ba7-49a2-be1b-5e4815ccb8f8
mysmatrix(TrivialDeltaLayer(), Diagonal([1, 1,1,1]), 1.0)

# ╔═╡ 09179a65-c8ed-448f-89fc-1aac43191bd0
struct TriviallyNonTrivialDeltaLayer <: DeltaScatterer end

# ╔═╡ 436afe39-699e-42b9-928f-868db66c2222
fft(Matrix(Diagonal([1, 1, 1])), 2)

# ╔═╡ ad8538e9-b49d-4390-9c4d-6fc4fe236fe1
Air(depth) = UniformLayer(depth, 1, 1)

# ╔═╡ 2e9cfa71-a698-4cd3-9785-0d5056b575b2
smatrix(Air(0), kxy, 1)[Block(2,1)]

# ╔═╡ 9efd0b2b-07d0-4949-b96f-88b01f78fa74
Matrix(I, 2,2 )

# ╔═╡ 1a7442fc-0944-4ea2-baf6-26a70dd9526b
function myTMsmatrix(layer::TriviallyNonTrivialDeltaLayer, kxy::NTuple{2, Frequencies}, ω::Real)
	kz = DeltaRCWA.get_kz(Air(1), kxy, ω)
	Nk = size(kz, 1)
	σₑ = zeros(Nk, Nk)
	σₘ = zeros(Nk, Nk)
	A = mortar(reshape([
		-I + σₑ * Diagonal(kz) / (2ω),
		-Diagonal(kz) / ω + σₘ/2,
		I - σₑ * Diagonal(kz) / (2ω),
		-Diagonal(kz) / ω + σₘ/2,
	], 2, 2))
	# return A
	B = mortar(reshape([
		I + σₑ * Diagonal(kz) / (2ω),
		-(Diagonal(kz) / ω + σₘ/2),
		-(I + σₑ * Diagonal(kz) / (2ω)),
		-(Diagonal(kz) / ω + σₘ/2),
	], 2, 2))
	return A\B
end

# ╔═╡ 1b63a590-2a4f-44fe-9ff8-31baac0eeb47
myTMsmatrix(TriviallyNonTrivialDeltaLayer(), (kx ,ky), 1.2)[Block(1,2)]

# ╔═╡ 21b38329-6cea-412e-837e-9038a44df9ea
kx

# ╔═╡ 15fea0cd-8d53-4413-bc3c-fd89cc73ec0e
DeltaRCWA.get_kz(Air(1), kxy, 1.0)

# ╔═╡ b4413619-b7bc-4909-89b3-21a22d2b0571
ifft(fft(Diagonal(1:4), 2), 1)

# ╔═╡ 92b9913c-de6d-4946-a82c-1e1568f3da5a


# ╔═╡ Cell order:
# ╠═a8253062-df3c-11eb-26e1-0dcff18bf886
# ╠═5c283082-518b-45c6-907a-e80945b1a4e4
# ╠═898fa516-591d-4831-8c65-f0e758922d15
# ╠═e944920d-7aa0-4fdc-ad47-4c49821eceeb
# ╠═ca120085-b234-48aa-9693-47171a799ea2
# ╠═7cc184e4-9f8d-4f2e-82cf-f5872b115102
# ╠═618bd4c7-a089-4d25-9491-16433c5479ae
# ╠═fc3390a9-caef-4602-b9fe-fe83cbe17cad
# ╠═f0580da9-52db-4c40-8706-b537870ab1dc
# ╠═4eb88559-8bbc-49f8-8d46-a4485892f57f
# ╠═159bdc89-6e78-435d-81e3-3cb2b77312e4
# ╠═153b9330-b0c4-4b84-9a26-aff06089b334
# ╠═02c7d3b2-f941-4973-8339-e5eb8bda0f19
# ╠═d8e567cb-90b8-4a36-902a-a980c684f728
# ╠═c2503f86-3083-4ebf-8a17-350bb1748e97
# ╠═25396236-c4ff-4161-8d1b-edb07174de2c
# ╠═dee4e053-cc7e-41ad-ba24-2a641131ce38
# ╠═9cf2ab2a-e9ab-4194-bdf3-5db2c4939a03
# ╠═4ff1f5e6-1cce-4da1-b2ef-99e491ea5719
# ╠═ade9e849-35ed-4943-9fd8-49f32d5fea61
# ╠═e89f082d-8cdb-4bc9-8d88-03e57c1fc11a
# ╠═91b692b2-bfff-4119-84ec-c3f6dc2e9120
# ╠═6e21e9a8-c7f0-474c-9730-2a0d4c768fb9
# ╠═63b477fd-38c5-4532-a752-916a5f08bee7
# ╠═53b9bb4a-96b9-40f5-85f5-42d6e0302703
# ╠═815d975c-0ad4-491d-bd2e-8f0053ffa29b
# ╠═fbd43517-8d23-4d53-8919-584a4b45817c
# ╠═4e78b99e-8084-4ed9-886d-f9bfd3190640
# ╠═801cf9ef-fe46-480a-9124-601fa03f62da
# ╠═fd1f4984-a0d0-498d-832f-1f884acf0265
# ╠═aaa1d158-59b3-4c72-9306-fb7bcef3a09e
# ╠═65ffffaa-6511-49df-b0ed-0353c52113d6
# ╠═3108c320-3f17-4e0b-858b-b9e25c11f02e
# ╠═484f4b6a-8ccb-4519-9d04-fb0b505fc13a
# ╠═f47a4fcb-fa16-48ab-850c-2a3de6bfeac7
# ╠═6d4dd4c0-8e0c-4825-aa64-3a82ced17cd3
# ╠═90d06519-9526-48c0-8560-116b7c112997
# ╠═9ba33324-8561-4413-8343-f9b52601d788
# ╠═0342e0f7-e654-41c2-9a56-b5e99d59ca27
# ╠═36a24a19-bbdc-466e-8905-da41493c5f8d
# ╠═7fdb6893-010e-4a6a-9636-d59e966e2fac
# ╠═3db65194-8647-41ff-8d30-8c79638b3753
# ╠═06a9e35c-0400-4875-9ea7-86927ff288fa
# ╠═315bffb4-ef0f-4dcf-8b0c-31fbd2f0a795
# ╠═7229e888-dc22-47d2-8dd5-db06bd5ae476
# ╠═9a7c15d1-6ac9-4e6b-851f-1ef502507254
# ╠═22081b08-f086-4a64-b349-b6938c1c4da1
# ╠═7a8ddcf3-40c4-4021-a197-8068d8cf22a2
# ╠═cc8c8ab2-d5b8-485a-a1e1-ce7e735bdf9c
# ╠═efe96bea-feee-413e-a03e-1254bc27e8d8
# ╠═ed5050f9-8f56-475c-9abd-d661a49aceb9
# ╠═1ef4e805-3dd9-4f39-a697-ddc0b4f04ad8
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
# ╠═2af29afd-37ea-4abf-b0c4-737cfd96e101
# ╠═a39dd1aa-868d-4819-877f-ccd4b856cb4c
# ╠═a022ef8b-03cc-42cc-ac9e-4da96e3edecb
# ╠═66a8554e-1b71-4776-a903-c9ee54f1d64b
# ╠═6045bf4a-e6ac-4411-bd0e-171ee51f3c97
# ╠═95cb80ae-4979-47e1-a798-11c8f33259f5
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
# ╠═992eb33a-3fcc-4d1b-bf48-3f33b8c2a1e4
# ╠═543c8328-ab82-435a-aca7-1f0ab573294d
# ╠═2e9cfa71-a698-4cd3-9785-0d5056b575b2
# ╠═8da39a74-6081-4bb1-8d6a-1a0d9385e189
# ╠═521f8552-1ba7-49a2-be1b-5e4815ccb8f8
# ╠═09179a65-c8ed-448f-89fc-1aac43191bd0
# ╠═436afe39-699e-42b9-928f-868db66c2222
# ╠═ad8538e9-b49d-4390-9c4d-6fc4fe236fe1
# ╠═9efd0b2b-07d0-4949-b96f-88b01f78fa74
# ╠═1a7442fc-0944-4ea2-baf6-26a70dd9526b
# ╠═1b63a590-2a4f-44fe-9ff8-31baac0eeb47
# ╠═21b38329-6cea-412e-837e-9038a44df9ea
# ╠═15fea0cd-8d53-4413-bc3c-fd89cc73ec0e
# ╠═b4413619-b7bc-4909-89b3-21a22d2b0571
# ╠═92b9913c-de6d-4946-a82c-1e1568f3da5a
