# Prototype Specification for scattering problems

(Drafted July 4, 2021)
(Comment July 22, 2021: interfaces still need to be defined and documented)

A scattering problem describes the effect of a process on a travelling medium,
such as when billiard balls collide or light waves pass through a lens.
Therefore the scattering process accepts inputs and returns outputs, usually
depending on some structure to the scattering interaction.
We already see that the medium which is travelling can be interacting with
other objects in the same medium (a fully-coupled scattering problem), or with
external structures which are typically assumed to scatter the medium without
being affected themselves (one-way coupling).

We can also make the distinction of a classical scattering problem (describing
probability distributions over the classical phase space of a particle)
from a quantum scattering problem (describing the amplitudes of wave-like
quantum fields in certain normal modes).

In practice, the mathematical structure of these different types of problems may
coincide in the sense that there is some medium (e.g. particles, fields)
interacting with an operator (e.g. potentials, collision operators, materials)
which removes the need for such granular distinctions.

Examples of a classical scattering process: billiard ball collisions,
Kolomogorov-type equations, Boltzmann kinetic theory, classical electromagnetism

Examples of a quantum scattering process: electron/neutron scattering

Scattering problems may also have other mathematical properties, such as
linearity (usually a result of the isotropy of space). I believe that the types
should emphasize the mathematical structure of the problem rather than the
physical details, so I am ok adding a LinearScatteringProblem type or a
IntegroDifferentialScatteringProblem type, but not a ClassicalScatteringProblem
type because the former may change the mathematics (and thus the methods) to
solve the problem whereas the latter is mostly for classification.

The preceding discussion mostly applies to the abstract types, each of which
may have its own set of generic solution methods, whereas the concrete types
may describe a physical type of problem, such as RCWAProblem.

A scattering problem involves a medium being scattered (such as billiard balls,
an electromagnetic field, a probability distribution) by a structure (such as an
object, a metalens, a potential) which is assumed to be stationary.


A scattering problem should be characterized by its dimensionality and its
symmetries. If each dimension is labeled by an axis, then an invariant axis
should simplify the problem by effectively eliminating one of the dimensions,
leading to a simpler solution method. Otherwise, periodicity along some axis
may provide a less expensive solution method than a problem without symmetries.

In every problem, there should be a domain over which to solve the problem

In summary, scattering problems have the following features:
- Geometry
    - dimensionality
    - symmetries
    - structures

The structure

# Prototype Specification for Layers

An AbstractScatteringLayer{N} is a type representing a linear scatterer in
N dimensions

Dimensionality conventions:
- 1D is along z, invariant in x, y
- 2D is along y, z, invariant in x or swap x and y
- 3D is along x, y, z
Any axis which is not invariant could still be periodic

A delta-function layer is measure zero along z and anything along x, y

All AbstractScatteringLayers should have a method to evaluate the scattering
matrix possibly by returning a dense matrix or in a matrix free formulation.
The layers should also have a well-defined impedance provided by a method, where
the impedance needs to be evaluated at some point. This provides the 
implmentation of the geometric model of the scattering layer which presumably
will be used internally to build the scattering matrix.

It is expected that each concrete type represents a layer with a parametrized
geometry, such as a Gaussian, or a unit cell of a periodic structure.
Then one has to write methods which give the impedance of this structure at all
points in this structure and methods to get the scattering matrix at a point.
Perhaps it is possible to write a generic method for the scattering matrix once
the impedance field is defined.
Instead of having a parametrized structure, one could also make a type that
specifies the impedances at grid points and interpolates the values elsewhere.

## Prototype names for methods
- SMatrix_? where ? takes values of: dense matvec
- 

## Prototype signatures for methods
SMatrix_? is going to evaluate the scattering matrix of a layer acting on a set
of modes, β, so its signature could be SMatrix_?{AbstractScatteringLayer, Modes}

## Subtypes

All 
Prototype names: DeltaScatteringLayer, ?ScatteringLayer where ? can be:
Thick, High, Deep, Material, Physical, Substantial, 3D, Volume
(other adjectives to describe a layer with measure > 0 volume)

UniformLayers are a concrete type that represent a homogenous slab of linear
dielectric and magnetic medium. They are translation invariant 
