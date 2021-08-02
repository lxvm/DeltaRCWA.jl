export RCWAProblem, solve

"""
    enforce_N_pol(N::Int64, ::P where P <: AbstractPolarization)

Check that a polarization is specified when there are symmetries which decouple
the TE and TM polarizations
"""
function enforce_N_pol(N::Int64, ::P) where P <: AbstractPolarization
    if 0 <= N <= 1
        @assert P <: UncoupledPolarization "1D, 2D photonic crystals uncouple TE, TM"
    elseif N == 2
        @assert P <: CoupledPolarization "3D photonic crystals couple TE, TM"
    else
        error("only 1D, 2D, 3D photonic crystals are possible")
    end
end

struct RCWAProblem{N}
    structure::RCWAScatterer{N, M} where M
    grid::NTuple{N, Tuple{Int64, Float64}}
    k⃗::NTuple{N, Frequencies}
    ω::T where T <: Real
    pol::P where P <: AbstractPolarization
    I₁::A where A <: AbstractArray{<: Number, N}
    I₂::B where B <: AbstractArray{<: Number, N}
    function RCWAProblem(
        structure::RCWAScatterer{N, M} where M,
        grid::NTuple{N, Tuple{Int64, Float64}},
        k⃗::NTuple{N, Frequencies},
        ω::T where T <: Real,
        pol::P where P <: AbstractPolarization,
        I₁::A where A <: AbstractArray{<: Number, N},
        I₂::B where B <: AbstractArray{<: Number, N},
    ) where N
        enforce_N_pol(N, pol)
        new{N}(structure, grid, k⃗, ω, pol, I₁, I₂)
    end
end

"""
    RCWAProblem(
        structure::RCWAScatterer{N, M} where M,
        grid::NTuple{N, Tuple{Int64, Float64}},
        ω::T where T <: Real,
        I₁::A where A <: AbstractArray{<: Number, N},
        I₂::B where B <: AbstractArray{<: Number, N},
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
    I₁::A where A <: AbstractArray{<: Number, N},
    I₂::B where B <: AbstractArray{<: Number, N},
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
    I₁::A where A <: AbstractArray{<: Number, N}
    I₂::B where B <: AbstractArray{<: Number, N}
    O₁::C where C <: AbstractArray{<: Number, N}
    O₂::D where D <: AbstractArray{<: Number, N}
    function RCWASolution(
        grid::NTuple{N, Tuple{Int64, Float64}},
        k⃗::NTuple{N, Frequencies},
        ω::T where T <: Real,
        pol::P where P <: AbstractPolarization,
        I₁::A where A <: AbstractArray{<: Number, N},
        I₂::B where B <: AbstractArray{<: Number, N},
        O₁::C where C <: AbstractArray{<: Number, N},
        O₂::D where D <: AbstractArray{<: Number, N},
    ) where N
        enforce_N_pol(N, pol)
        new{N}(grid, k⃗, ω, pol, I₁, I₂, O₁, O₂)
    end
end

function solve(P::RCWAProblem)::RCWASolution
    @assert all(size(c) == Tuple(e[1] for e in P.dims) for c in (P.C₁, P.C₂))
    smatrix(P.structure, P.modes) * to_blockvector(P.incident)
    # TODO: test these properties of the solution:
    # Nondimensionalize the problem in terms of the lattice vectors
    # fftfreq for each spatial grid
    # construct the inital mode as a blockvector
end

function to_blockvector(input::NTuple{2, NTuple{M, ComplexF64}}) where M
    BlockVector(collect([e for block in input for e in block]), [M, M])
end