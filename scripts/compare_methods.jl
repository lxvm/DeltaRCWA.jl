using Pkg
Pkg.activate(ENV["HOME"] * "/.julia/dev/DeltaRCWA/")

using DeltaRCWA

### Carlos' solver

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

# u_p = ComplexF64[n == 0 ? 1 : 0 for n in nvec]
u_p = fft(uInc(xvec, 0))/length(nvec)
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

E1(x,z) = sum([u_out_p[n+1]*exp(-1im*kₓ[n+1]*x)*exp(1im*β[n+1]*z) for n in nvec])+sum([u_p[n+1]*exp(-1im*kₓ[n+1]*x)*exp(1im*β[n+1]*z) for n in nvec])
E2(x,z) = sum([u_out_n[n+1]*exp(-1im*kₓ[n+1]*x)*exp(-1im*β[n+1]*z) for n in nvec])

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

### DeltaRCWA solver
struct ComplexExpSheet{T <: Real} <: RCWASheet{1}
    θ::T
    θᵗ::T
    k::T
    d::T
    L::T
end

function σₑˣˣ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
	2 ./ M₀(x⃗...)
	# ComplexF64[1 - im*sin(sheet.k*sheet.d*x[1])/(1+cos(sheet.k*sheet.d*x[1])) for x in Iterators.product(x⃗...)]
end

function σₘʸʸ(sheet::ComplexExpSheet{T}, x⃗::Tuple{StepRangeLen}) where T <: Real
	2N₀(x⃗...)
end

ω = k
sheet = ComplexExpSheet(θ, θᵗ, k, d, L)
mode_basis = PlanewaveModes(ω, ((length(nvec), sheet.L), ), Vacuum())
pol = TM()
prob = DeltaRCWAProblem(sheet, mode_basis, pol, u_p, u_n)
sol = solve(prob)