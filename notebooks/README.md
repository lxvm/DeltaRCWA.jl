# Notebooks

`DeltaRCWA.ipynb` is a Jupyter notebook exploring the scattering matrix of a sheet-transition condition.

`star_test.jl` is a Pluto notebook explaining `DeltaRCWA.ipynb` and providing some comparisons to boundary integral techniques for the same

`interface.jl` describes how to use the interface provided by `DeltaRCWA` and compares the results to a direct numerical calculation.

`derivations.jl` uses `Symbolics.jl` to derive the scattering matrix from Generalized Sheet Transition Conditions.

`inversion.jl` demonstrates use of iterative solvers, compares their accuracy
at inverting blocks of a matrix, and benchmarks different methods.

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

Now the environment should be activated and starting one of the notebooks can be done with `using Pluto; Pluto.run()` or `using IJulia; notebook()` depending on the type of notebook.

Instructions
- Install [Julia](https://julialang.org/)
- Clone [`lxvm/DeltaRCWA.jl`](https://github.com/lxvm/DeltaRCWA.jl) into the folder `~/.julia/dev/DeltaRCWA` with `git`
- `cd` to `~/.julia/dev/DeltaRCWA/notebooks` folder and start a `julia` session
- Activate the package manager and notebook environment by typing `]activate .` in the REPL prompt
- Install DeltaRCWA.jl in this environment at the package prompt with `dev DeltaRCWA`
- Finish installing the environment at the package prompt with `instantiate`