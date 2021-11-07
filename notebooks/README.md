# Notebooks

`3Dproblem.jl` verifies the 3D scattering matrix against the 2D one.

`interface.jl` describes how to use the interface provided by `DeltaRCWA`.

`GSTC.jl` uses `Symbolics.jl` to form scattering matrix from GSTCs.

`inversion.jl` demonstrates use of iterative solvers, compares their accuracy
at inverting blocks of a matrix, and benchmarks different methods.

`Fresnel.jl` shows how GSTCs recover Fresnel scattering.

## Setting up this environment

TL;DR
```
$ git clone https://github.com/lxvm/DeltaRCWA.jl.git ~/.julia/dev/DeltaRCWA
$ cd ~/.julia/dev/DeltaRCWA/notebooks
$ julia
julia>]
pkg> activate .
pkg> dev DeltaRCWA
pkg> instantiate
```

Now the environment should be activated and starting one of the notebooks can be
done with `using Pluto; Pluto.run()`.

Instructions
- Install [Julia](https://julialang.org/)
- Clone [`lxvm/DeltaRCWA.jl`](https://github.com/lxvm/DeltaRCWA.jl) into the folder `~/.julia/dev/DeltaRCWA` with `git`
- `cd` to `~/.julia/dev/DeltaRCWA/notebooks` folder and start a `julia` session
- Activate the package manager and notebook environment by typing `]activate .` in the REPL prompt
- Install DeltaRCWA.jl in this environment at the package prompt with `dev DeltaRCWA`
- Finish installing the environment at the package prompt with `instantiate`