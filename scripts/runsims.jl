using JLD2

include("Carlos_method.jl")
include("Luke_method.jl")
include("Lorenzo_method.jl")

include("analyses/show.jl")
include("analyses/interpolate.jl")
include("analyses/errors.jl")

## Define metasurface examples

function grating(; L=1.0, f=0.1)
    L₀ = L # unit cell width
    η(x) = zero(ComplexF64)
    μ(x) = abs(x - L₀/2) < L₀*f/2 ? im*one(ComplexF64) : zero(ComplexF64) # im in post, 0 outside
    (name="grating", L=L, η=η, μ=μ)
end

function deflector(; L=1.0, θⁱ=π/2, θᵗ=π/3)
    # these parameters are taken from https://doi.org/10.1364/OE.26.030202 
    L₀ = L # unit cell width
    d = cos(θᵗ)-cos(θⁱ) # design parameter
    k₀ = 2π/(L₀*abs(d)) # wavenumber that attains the metasurface design goals
    η(x) = sin(abs(θⁱ))*(1+exp(im*k₀*d*x))
    μ(x) = sin(abs(θⁱ))*(1-exp(im*k₀*d*x))
    (name="deflector", L=L, η=η, μ=μ)
end

function bump(; L=1.0, A=1.0, σ=0.1)
    L₀ = L # unit cell width
    σ₀ = σ/L₀ # normalized standard deviation
    η(x) = zero(ComplexF64)
    μ(x) = (A/(σ₀*sqrt(2π)))*exp(-0.5((x-L₀/2)/σ₀)^2)
    (name="bump", L=L, η=η, μ=μ)
end

function imbump(; L=1.0, A=1.0, σ=0.1)
    L₀ = L # unit cell width
    σ₀ = σ/L₀ # normalized standard deviation
    η(x) = zero(ComplexF64)
    μ(x) = im*(A/(σ₀*sqrt(2π)))*exp(-0.5((x-L₀/2)/σ₀)^2)
    (name="imbump", L=L, η=η, μ=μ)
end

## Define simulation and analysis sweeps

# sim should be a NamedTuple with fields "method", "param", "sheet", and
# optionally "sol"

# param should be a NamedTuple with fields "n", "k", "i" corresponding to the
# number of grid points in unit cell, a multiplier for the wavenumber of the
# metasurface period, and the index of the incident Bloch mode

function name_sim(sim)
    (;method, sheet, param) = sim
    string(sheet.name, "/L", sheet.L, "/", param.k, "k₀/i", param.i, "/n", param.n, "/", name_method(method))
end

name_method(method) = split(string(typeof(method).name.name), "_")[2]

function run_sim(method, param, sheet)
    sim = (method=method, param=param, sheet=sheet)
    k₀ = 2π / sheet.L
    @info string("Running ", name_sim(sim))
    @time sol = method(sheet.η, sheet.μ, sheet.L, param.n, k₀ * param.k, param.i)
    merge(sim, (sol=sol, ))
end

function gen_sweep(methods, params, sheets)
    (run_sim(m, p, s) for (m, p, s) in Iterators.product(methods, params, sheets))
end

function gen_names(methods, params, sheets)
    (name_sim((method=m, param=p, sheet=s)) for (m, p, s) in Iterators.product(methods, params, sheets))
end

function gen_params(e, k, i; b=2)
    ((n=p[1], k=p[2], i=p[3]) for p in Iterators.product(b .^ e, k, i))
end

function save_sim(file, name, method, param, sheet)
    if haskey(file, name)
        @info string("already found ", name)
        sim = read_sim(file, name)
    else
        sim = run_sim(method, param, sheet)
        file[name] = sim.sol
    end
    return sim
end

read_sim(file, name) = merge(recover_sim(name), (sol=file[name],))

function recover_sim(name)
    parts = split(name, "/")
    method = eval(Symbol("compute_" * parts[end] * "_method"))
    sheet = eval(Symbol(parts[1]))(; L=parse(Float64, parts[2][2:end]))
    param = (n=parse(Int, parts[5][2:end]), k=parse(Float64, parts[3][1:(length(parts[3])-2)]), i=parse(Int, parts[4][2:end]))
    (method=method, param=param, sheet=sheet)
end

function save_sweep(f, m, p, s; filename="results.jld2")
    jldopen(filename, "a+") do file
        f(save_sim(file, x[1], x[2]...) for x in zip(gen_names(m, p, s), Iterators.product(m, p, s)))
    end
end

function read_sweep(f, methods, params, sheets; filename="results.jld2")
    jldopen(filename, "a+") do file
        f(read_sim(file, name) for name in gen_names(methods, params, sheets))
    end
end

function make_param_plot(sim; prefix="analyses/")
    path = prefix * join(split(name_sim(sim), "/")[1:2], "/") * "/param.png"
    mkpath(dirname(path))
    savefig(show_params(sim), path)
end

function make_field_plot(sim; prefix="analyses/")
    path = prefix * name_sim(sim) * "field.png"
    mkpath(dirname(path))
    savefig(show_fields(to_position_basis(sim)), path)
end

function make_convergence_plot(x; b=2, prefix="analyses/")
    L2errors = Vector{NamedTuple{(:I₁, :O₁, :I₂, :O₂), NTuple{4, Float64}}}(undef, Int(length(x)//2))
    exponents = Vector{Float64}(undef, Int(length(x)//2))
    for (i, e) in enumerate(Iterators.partition(x, 2))
        sims = to_position_basis.(collect(e))
        L2errors[i] = errorL2(get_errors(sims[1], sims[2]))
        exponents[i] = log2(sims[1].param.n)
    end
    path = prefix * join(split(name_sim(first(x)), "/")[1:4], "/") * "/convergence.png"
    mkpath(dirname(path))
    savefig(show_errorL2(exponents, L2errors), path)
end
function make_convergence_plots(x)
    dims = size(x) # 1st dim is methods, then n, k, i, examples
    dims[1] == 2 || error("unexpected number of methods: compare exactly two")
    ndims(x) == 5 || error("unexpected ndims: supply all parameters in a vector")
    make_convergence_plot.(Iterators.partition(x, Int(length(x)//prod(dims[3:5]))))
end

function make_distance_plot(x; prefix="analyses/")
    sims = to_position_basis.(collect(x))
    path = prefix * dirname(name_sim(first(x))) * "/errordistance.png"
    mkpath(dirname(path))
    savefig(show_errorbydistance(get_errors(sims[1], sims[2])), path)
end
function make_distance_plots(x)
    dims = size(x) # 1st dim is methods, then n, k, i, examples
    dims[1] == 2 || error("unexpected number of methods: compare exactly two")
    ndims(x) == 5 || error("unexpected ndims: supply all parameters in a vector")
    make_distance_plot.(Iterators.partition(x, Int(length(x)//prod(dims[2:5]))))
end

function make_error_plot(x; prefix="analyses/")
    sims = to_position_basis.(collect(x))
    path = prefix * dirname(name_sim(first(x))) * "/errorfield.png"
    mkpath(dirname(path))
    savefig(show_fields(get_errors(sims[1], sims[2]),), path)
end
function make_error_plots(x)
    dims = size(x) # 1st dim is methods, then n, k, i, examples
    dims[1] == 2 || error("unexpected number of methods: compare exactly two")
    ndims(x) == 5 || error("unexpected ndims: supply all parameters in a vector")
    make_error_plot.(Iterators.partition(x, Int(length(x)//prod(dims[2:5]))))
end

# typical commands
# save_sweep(x -> make_field_plot.(x), [compute_BIE_method, compute_DeltaRCWA_method], gen_params(5:8, [2.2, 1.8], 1), [imbump(), ])
# read_sweep(make_convergence_plots, [compute_BIE_method, compute_DeltaRCWA_method], gen_params(5:8, [2.2, 1.8], [1]), [imbump(), ])
# read_sweep(make_distance_plots, [compute_BIE_method, compute_DeltaRCWA_method], gen_params(5:8, [2.2, 1.8], [1]), [imbump(), ])
# read_sweep(make_error_plots, [compute_BIE_method, compute_DeltaRCWA_method], gen_params(5:8, [2.2, 1.8], [1]), [imbump(), ])
# read_sweep(x -> make_param_plot(first(x)), [compute_BIE_method, compute_DeltaRCWA_method], gen_params(5:8, [2.2, 1.8], [1]), [imbump(), ])