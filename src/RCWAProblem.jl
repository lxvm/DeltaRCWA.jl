# TODO: finish creating the flexible interface for arbitrary RCWA problems
# export RCWAProblem, solve

struct RCWAProblem{N}
    structure::RCWAScatterer{N, M} where M
    grid::NTuple{N, Tuple{Int64, Float64}}
    k⃗::NTuple{N, Frequencies}
    ω::Float64
    pol::AbstractPolarization
    I₁::WeightedModes{N}
    I₂::WeightedModes{N}
    function RCWAProblem(
        structure::RCWAScatterer{N, M} where M,
        grid::NTuple{N, Tuple{Int64, Float64}},
        k⃗::NTuple{N, Frequencies},
        ω::T where T <: Real,
        pol::P where P <: AbstractPolarization,
        I₁::WeightedModes{N},
        I₂::WeightedModes{N},
    ) where N
        enforce_N_pol(N, pol)
        new{N}(structure, grid, k⃗, convert(Float64, ω), pol, I₁, I₂)
    end
end

"""
    RCWAProblem(
        structure::RCWAScatterer{N, M} where M,
        grid::NTuple{N, Tuple{Int64, Float64}},
        ω::T where T <: Real,
        I₁::WeightedModes{N},
        I₂::WeightedModes{N},
        pol::P where P <: AbstractPolarization,
    ) where N

Stores information needed to solve a scattering problem: the structure to
scatter from, the dimensionality of the problem, the incident modes to
scatter and their frequency

`grid` is tuple of Tuple{Int64, Float64} describing each periodic axis
by a number of grid points (N), the lattice constant (a).
This will determine the spatial grid in each direction == 0:LC/N:LC .- origin.

`ω` is the frequency of the incident wave

`I₁, I₂` are of type `AbstractArray` and contain the mode amplitudes of the
incident waves.

`pol` is an concrete subtype of `AbstractPolarization`
"""
function RCWAProblem(
    structure::RCWAScatterer{N, M} where M,
    grid::NTuple{N, Tuple{Int64, Float64}},
    ω::T where T <: Real,
    pol::P where P <: AbstractPolarization,
    I₁::WeightedModes{N},
    I₂::WeightedModes{N},
) where N
    k⃗ = Tuple(fftfreq(e[1], e[1] / e[2]) for e in grid)
    RCWAProblem(structure, grid, k⃗, ω, pol, I₁, I₂)
end

"""
Store the results of a RCWA solve, namely the incident and scattered mode
amplitudes and information about the modes + spatial grid in order to
reconstruct the total fields.
"""
struct RCWASolution{N}
    grid::NTuple{N, Tuple{Int64, Float64}}
    k⃗::NTuple{N, Frequencies}
    ω::T where T <: Real
    pol::P where P <: AbstractPolarization
    I₁::WeightedModes{N}
    I₂::WeightedModes{N}
    O₁::WeightedModes{N}
    O₂::WeightedModes{N}
    function RCWASolution(
        grid::NTuple{N, Tuple{Int64, Float64}},
        k⃗::NTuple{N, Frequencies},
        ω::T where T <: Real,
        pol::P where P <: AbstractPolarization,
        I₁::WeightedModes{N},
        I₂::WeightedModes{N},
        O₁::WeightedModes{N},
        O₂::WeightedModes{N},
    ) where N
        enforce_N_pol(N, pol)
        new{N}(grid, k⃗, ω, pol, I₁, I₂, O₁, O₂)
    end
end

function solve(P::RCWAProblem)::RCWASolution
    @assert all(size(c) == Tuple(e[1] for e in P.dims) for c in (P.C₁, P.C₂))
    smatrix(P.structure) * to_blockvector(P.incident)
    # TODO: test these properties of the solution:
    # Nondimensionalize the problem in terms of the lattice vectors
    # fftfreq for each spatial grid
    # construct the inital mode as a blockvector
end

function to_blockvector(input::NTuple{2, NTuple{M, ComplexF64}}) where M
    BlockVector(collect([e for block in input for e in block]), [M, M])
end