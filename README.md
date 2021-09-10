# DeltaRCWA.jl

An approximate frequency-domain solver for electromagnetic fields using
[RCWA](https://en.wikipedia.org/wiki/Rigorous_coupled-wave_analysis).

## Contents

- `notebooks`: interactive Pluto and Jupyter notebooks showcasing the uses of
`DeltaRCWA` as well as doing exploratory tests/debugging which will hopefully
be ready to include in `DeltaRCWA`
- `src`: the source code of `DeltaRCWA`
- `scripts`: contains Julia scripts to do convergence and validation tests
for different solution methods
- `test`: contains unit tests of `DeltaRCWA`'s star product methods

The root folder and `notebooks`, `scripts`, and `test` each contain a separate
Julia environment. See `notebooks/README.md` for an example on how to set one up.

## Documentation

Currently the documentation is informal at its best.
An attempt has been made to provide docstrings for most types and functions.
To see an example usage, open the Pluto notebook in `notebooks/interface.jl`.

## Tests

This is todo and not to be done until more of the package is working.