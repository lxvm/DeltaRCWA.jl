script_dir = ENV["HOME"] * "/.julia/dev/DeltaRCWA/scripts/"
import Pkg
Pkg.activate(script_dir)

using FFTW
using LinearAlgebra
using Plots
using JLD2

include("Carlos_method.jl")
include("Luke_method.jl")
include("Lorenzo_method.jl")

η(x) = 0
μ(x) = 0
L=1.0
N=100
k=11.0
i=1

println("starting BIE")
BIE = compute_BIE_method(η, μ, L, N, k, i)
println("starting RCWA")
RCWA = compute_RCWA_method(η, μ, L, N, k, i)
println("starting DeltaRCWA")
ΔRCWA = compute_DeltaRCWA_method(η, μ, L, N, k, i)