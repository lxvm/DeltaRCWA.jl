### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 83d9578c-3960-486c-ad3d-b4449f7bae6f
md"
# Fresnel scattering at a dielectric interface

Note that our GSTCs simplify to continuity of the tangential fields when
the surface conductivities vanish.
Therefore, our scattering matrix should recover the Fresnel coefficients
when we allow the medium on each side of the interface to differ.
"

# ╔═╡ 325b8359-ebac-434b-bba1-31aba05c2c39
md"
## Higher order corrections for a metasurface at a dielectric interface

The GSTCs we use are published by Epstein (DOI: 10.1109/TAP.2014.2354419)

These are based on ones published by Kuester (DOI: 10.1109/TAP.2003.817560), which
make the following assumptions.
'In this paper, we will derive what we believe to be the most flexible form of a
generalized sheet transition condition (GSTC) applicable to metafilms located in a
homogeneous medium. The case of a metafilm on a dielectric interface will be
considered in a separate paper.'

Namely, we shouldn't expect to recover Fresnel scattering by simply using the same
GSTC with different impedance parameters for each film.
Kuester published a follow-up article (DOI: 10.1109/APS.2010.5562250)
with a first order correction to the previous one which should generalize to the
case where the metafilm has different media on either side.
"

# ╔═╡ Cell order:
# ╠═83d9578c-3960-486c-ad3d-b4449f7bae6f
# ╠═325b8359-ebac-434b-bba1-31aba05c2c39
