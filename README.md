# DeltaRCWA.jl

An approximate frequency-domain solver for electromagnetic fields using
[RCWA](https://en.wikipedia.org/wiki/Rigorous_coupled-wave_analysis).

## Contents

- `notebooks`: interactive Pluto and Jupyter notebooks showcasing the functionality
of `DeltaRCWA` as well as doing exploratory tests/debugging which will hopefully
be ready to include in `DeltaRCWA`
- `src`: the source code of `DeltaRCWA`
- `scripts`: contains Julia scripts to do integration tests of different solution
methods
- `test`: contains unit tests of `DeltaRCWA`'s star product methods

The root folder and `notebooks`, `scripts`, and `test` each contain a separate
Julia environment. See `notebooks/README.md` for an example on how to set one up.