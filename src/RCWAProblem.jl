using Base: Real
"""
    RCWAProblem(
        structure::AbstractScatterer,
        a::NTuple{N, Float64} where N,
        ω::Float64,
        incidentmodes::NTuple{N, NTuple{2, NTuple{M, ComplexF64} where M}},
    )

Stores information needed to solve a scattering problem: the structure to
scatter from, the dimensionality of the problem, the incident modes to
scatter and their frequency


Can be used to solve a scattering problem for Maxwell's equations in
N + 1 dimensions where 0 <= N <= 2. The additional dimension is the propagation
axis where the fields are calculated, and the N dimensions are periodic.
Any remaining dimensions are assumed to be invariant

`a` is tuple of the lattice constants for each periodic axis
`ω` is the frequency of the incident wave
`incidentmodes` is a Tuple containing the incident mode coefficients
"""
struct RCWAProblem{N} <: PeriodicScatteringProblem where N
    structure::AbstractScatterer
    a::NTuple{N, Float64}
    ω::Float64
    incidentmodes::NTuple{N, NTuple{2, NTuple{M, ComplexF64} where M}}
end

function solve(problem::RCWAProblem{N}) where N
    error("Not implemented")
    return smatrix(
        problem.structure,
        (),
        problem.ω) * blockvector_product(problem.incidentmodes)
    # TODO
    # Nondimensionalize the problem in terms of the lattice vectors
    # fftfreq for each spatial grid
    # construct the inital mode as a blockvector
end

function blockvector_product(in::NTuple{N, NTuple{M, ComplexF64} where M} where N)
    #TODO
end