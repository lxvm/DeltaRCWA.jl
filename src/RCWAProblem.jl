"""
    RCWAProblem(
        structure::AbstractScatterer,
        a::NTuple{N, Float64} where N,
        ω::Float64,
        incidentmodes::NTuple{2, NTuple{M, ComplexF64} where M}
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
struct RCWAProblem{N, P} <: PeriodicScatteringProblem
    structure::DeltaScatterer
    modes::PlanewaveModes{N, P, T} where T <: AbstractUniformMedium
end

function solve(problem::RCWAProblem{N, P}) where {N, P}
    return smatrix(
        problem.structure,
        Tuple(fftfreq(length(problem.modes[i][1]), 1/a[i]) for i in eachindex(a)),
        problem.ω
    ) * to_blockvector(problem.modes)
    # TODO: test these properties of the solution:
    # Nondimensionalize the problem in terms of the lattice vectors
    # fftfreq for each spatial grid
    # construct the inital mode as a blockvector
end

function to_blockvector(input::NTuple{2, NTuple{M, ComplexF64}}) where M
    return BlockVector(
        collect([e for block in input for e in block]),
        [M, M]
    )
end