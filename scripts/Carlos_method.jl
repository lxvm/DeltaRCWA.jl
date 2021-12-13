nbdir = ENV["HOME"] * "/.julia/dev/DeltaRCWA/notebooks/"
include("$(nbdir)NystromMethodQP.jl")
import FFTW: fftfreq
import LinearAlgebra: I, Diagonal

"""
    compute_BIE_method(η::Function, μ::Function, L::Float64, N::Int, k::Float64, i::Int)

Compute the fields using Carlos' Boundary Integral Equation (BIE) method.
The output fields are on a grid compatible with the RCWA method.

Arguments:
η:: A periodic function returning the metasurface magnetic conductance
μ:: A periodic function returning the metasurface electric resistance
L:: the length of the unit cell (some multiple of the period of η, μ)
N:: the number of discretization points along the unit cell
k:: the wavenumber/frequency of the incident planewave (assumed TM polarization)
i:: the index of the incident planewave in the Fourier basis

Keyword arguments:
method:: the Nystrom method to be used (choose default "MK" or "Alpert)
NC:: the number of cells in a windowed sum for Green's function (default 5)
NP:: the number of points in BIE discretization per wavelength (default 20)
s:: a scale factor to stretch the y coordinate of the output grid

Returns:
NamedTuple with fields:
I₁:: an array of the incident fields in the position basis for y<0
O₁:: an array of the outgoing fields in the position basis for y<0
I₂:: an array of the incident fields in the position basis for y>0
O₂:: an array of the outgoing fields in the position basis for y>0
"""
function compute_BIE_method(
    η::Function, μ::Function, L::Float64, N::Int, k::Float64, i::Int;
    method="MK", NC=5, NP=20, s=1
)
    dx = L/N
    k⃗ = 2π*fftfreq(N, 1/dx)
    β⃗ = @. sqrt(Complex(k^2 - k⃗^2))
    α = k⃗[i]
    β = β⃗[i]
    uInc(x,y)= exp(im*α*x+im*β*y)  # incident planewave
    # Defining geometry
    # Vertices that define the curve Γ 
    ver = [
        L/2   0;   #1 
        -L/2   0   #2
    ]
    # In this case, Γ is just a line segment from (L/2,0) to (-L/2,0). 
    # This matrix defines the orientation of the curve
            #starting node  #ending node  #domain on the right #domain on the left  
    info = [1               2             1                    0]
    # The parameter M defines the number of boundary points to used in the discretization of the BIE
    λ = 2π/k
    M  = Int(round(max(NP,NP/λ))) # ~20 pts per unit length
    # Γ contains all the information of the discretized curve that is needed for the numerical solution of the problem
    Γ  = SkeletonCctor(ver,info,[M])

    # Solving the boundary integral equations
    # method::Selection of the Nystrom method to be used
    # NC::Number of unit cells considered in the windowed sum (A = NC*L)
    mat = matricesSkeleton(1,Γ,k,L,α,NC,method) # Construction of the Nystrom matrices
        
    x1 = Γ.parts[1].x[:,1] # x-coordinate of the discretization points on Γ

    f1(x) = -2*im*k*η(x)*exp(im*α*x)
    A⁺ = -0.5*I+((im*k)*Diagonal(η.(x1)))*mat.S  # System matrix for Ω⁺
    b⁺ = f1.(x1) # boundary data for Ω⁺
    ψ⁺ = A⁺\b⁺

    f2(x) = -2*im*β*exp(im*α*x)
    A⁻ = -0.5*I+((im*k)*Diagonal(μ.(x1)))*mat.S # System matrix for Ω⁻
    b⁻ =  f2.(x1) # boundary data for Ω⁻
    ψ⁻ = A⁻\b⁻

    φ⁺ = 0.5*(ψ⁺+ψ⁻)
    φ⁻ = 0.5*(ψ⁺-ψ⁻)

    x = range(0, step=dx, length=N)
    y = range(-s*L, s*L, length=2N)
    y₁ = y[y .< 0] # y points in domain below sheet
    y₂ = y[y .> 0] # y points in domain above sheet

    ver⁺ = [L 0;-L 0;-L 2s*L;L 2s*L]
    Ω⁺ = isin(x,y,ver⁺)

    ver⁻ = [L 0;-L 0;-L -2s*L;L -2s*L]
    Ω⁻ = isin(x,y,ver⁻)

    UInc = reshape(uInc.(Ω⁺.pts[:,1],Ω⁺.pts[:,2]),Ω⁺.Ny,Ω⁺.Nx)

    pot⁺ = potentialsSkeleton(Ω⁺.pts,Γ,k,L,α,NC);
    U⁺scat = reshape(pot⁺.SL*φ⁺,Ω⁺.Ny,Ω⁺.Nx)
    U⁺   =  Ω⁺.In .* (U⁺scat .+ UInc)

    pot⁻= potentialsSkeleton(Ω⁻.pts,Γ,k,L,α,NC);
    U⁻scat = reshape(pot⁻.SL*φ⁻,Ω⁻.Ny,Ω⁻.Nx)
    U⁻   =  Ω⁻.In .* (U⁻scat .+ UInc)
    
    Ω⁺mask = isfinite.(Ω⁺.In)
    I₂ = reshape(UInc[Ω⁺mask] .- UInc[Ω⁺mask], length.((y₂, x)))
    O₂ = reshape(U⁺[Ω⁺mask], length.((y₂, x)))

    Ω⁻mask = isfinite.(Ω⁻.In)
    I₁ = reshape(UInc[Ω⁻mask], length.((y₁, x)))
    O₁ = reshape(U⁻[Ω⁻mask] .- UInc[Ω⁻mask], length.((y₁, x)))
    (x=x, y₁=y₁, y₂=y₂, I₁=I₁, O₁=O₁, I₂=I₂, O₂=O₂)
end