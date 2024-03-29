### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ c24044bc-2fc3-4e34-b99c-7dc4ec6a7515
import Pkg

# ╔═╡ 42da96ce-8122-4aa0-99d6-635f68db6e66
Pkg.activate(".")

# ╔═╡ 88e46392-1fbc-4079-8207-dd5b3ac01d7d
using DeltaRCWA

# ╔═╡ 16316725-c608-49a5-b90c-8caff973683a
using Plots

# ╔═╡ 83d9578c-3960-486c-ad3d-b4449f7bae6f
md"
# Fresnel scattering at a dielectric interface

Note that our GSTCs simplify to continuity of the tangential fields when
the surface conductivities vanish.
Therefore, our scattering matrix should recover the Fresnel coefficients
when we allow the medium on each side of the interface to differ.

The analytical solution for this problem is known and I will just cite it [here](https://homerreid.github.io/scuff-em-documentation/tests/FresnelScattering/FresnelScattering/)
``\newcommand{\bm}{\boldsymbol}``
"

# ╔═╡ 6cbc5f34-fe14-45e8-9fee-9787d76ad0fe
struct EmptySheet <: Sheet end

# ╔═╡ fb8df44c-a354-42de-9e99-df64dd4c5727
DeltaRCWA.ElectricResponseStyle(::Type{EmptySheet}) = Admittance()

# ╔═╡ f3d36c35-c09b-483f-bb41-046b6ecd0773
const Water = UniformMedium{1.7, 1.0}()

# ╔═╡ e899dbc0-03bb-4a13-9c09-bfad33e5c1d3
const stack = SheetStack(EmptySheet(), (Vacuum, Water))

# ╔═╡ 43d29086-510d-4870-94b7-2ff5ff34324b
n = 2 # number of grid points

# ╔═╡ 5377ae36-32fc-41c5-bd03-c737422b17f4
prob = DeltaRCWAProblem(stack, ((n, 1.0), ), 10.0, [i==2 ? 1 : 0 for i in 1:n], zeros(n))

# ╔═╡ 65c2e966-8cb0-4d76-96db-18f9045695a3
sol = solve(prob)

# ╔═╡ bf1213cd-2ec1-4d1e-8ebc-5be5af26621c
plot(sol; combine=true)

# ╔═╡ ba0c4060-f9d4-43b6-bc31-b48dcd81dd91
md"""
``k_z = \sqrt{\omega^2 \epsilon \mu - k_x^2 - k_y^2}``
"""

# ╔═╡ 687627f9-6938-4d7a-9395-9f1783b8a449
md"""
## Scattering matrix in 2d

We restrict ourselves to the decoupled TE/TM mode case (2d with no off-diagonal impedance matrix elements) for simplicity (we can make the 3d case work too).
Starting from the GSTCs for the TM modes (when ``k_y=0`` and the TE mode zero, ``E_y=H_x=0``)
```math
\begin{align}
-(H_y^{(2,+)} + H_y^{(2,-)} - H_y^{(1,+)} - H_y^{(1,-)})
	&= \sigma^e_{xx} (E_x^{(2,+)} + E_x^{(2,-)} + E_x^{(1,+)} + E_x^{(1,-)})/2
\\
(E_x^{(2,+)} + E_x^{(2,-)} - E_x^{(1,+)} - E_x^{(1,-)})
	&= -\sigma^m_{yy} (H_y^{(2,+)} + H_y^{(2,-)} + H_y^{(1,+)} + H_y^{(1,-)})/2
\end{align}
```
and the corresponding planewave relations (indexed by medium ``i``)
```math
\begin{align}
% -\omega^2 \epsilon^i \mu^i H_x^i &= \omega \epsilon^i k_z^i E_y^i
% \\
(k_x^2 -\omega^2 \epsilon^i \mu^i) H_y^i =  -(k_z^i)^2 H_y^i &= -\omega \epsilon^i k_z^i E_x^i
\\
-\omega^2 \epsilon^i \mu^i E_x^i &= -\omega \mu^i k_z^i H_y^i
% \\
% (k_x^2 -\omega^2 \epsilon^i \mu^i) E_y^i &= \omega \mu^i k_z^i H_x^i
\end{align}
```
which both are the same and let us eliminate the electric field using
```math
\begin{align}
E_x^i &= \frac{k_z^i}{\omega \epsilon^i} H_y^i
\end{align}
```
Taking into account the direction of propagation means
```math
\begin{align}
E_x^{(i,\pm)} &= \pm\frac{k_z^i}{\omega \epsilon^i} H_y^{(i,\pm)}
\end{align}
```
Note that this is the point where the 2d case is advantageous because the decoupled modes make the inversion in this case just the inverse of a diagonal matrix.
Hence the GSTCs are
```math
\begin{align}
&- \left( H_y^{(2,+)} + H_y^{(2,-)} - H_y^{(1,+)} - H_y^{(1,-)} \right)
\\ &\qquad = \sigma^e_{xx} \left( \frac{k_z^2}{2\omega\epsilon^2} \left( H_y^{(2,+)} - H_y^{(2,-)} \right) + \frac{k_z^1}{2\omega\epsilon^1} \left( H_y^{(1,+)} - H_y^{(1,-)} \right) \right)
\\
&\frac{2k_z^2}{\omega\epsilon^2} \left( H_y^{(2,+)} - H_y^{(2,-)} \right) - \frac{2k_z^1}{\omega\epsilon^1} \left( H_y^{(1,+)} - H_y^{(1,-)} \right)
\\ &\qquad = -\sigma^m_{yy} \left( H_y^{(2,+)} + H_y^{(2,-)} + H_y^{(1,+)} + H_y^{(1,-)} \right)
\end{align}
```
which gives the following scattering matrix
```math
\begin{align}
&\begin{pmatrix}
I + \sigma^e_{xx} \frac{k_z^1}{2\omega\epsilon^1}
&
-I - \sigma^e_{xx} \frac{k_z^2}{2\omega\epsilon^2}
\\
\sigma^m_{yy} + \frac{2k_z^1}{\omega\epsilon^1}
&
\sigma^m_{yy} + \frac{2k_z^2}{\omega\epsilon^2}
\end{pmatrix}
\begin{pmatrix}
H_y^{(1,-)}
\\
H_y^{(2,+)}
\end{pmatrix}
\\ =
&\begin{pmatrix}
-I + \sigma^e_{xx} \frac{k_z^1}{2\omega\epsilon^1}
&
I - \sigma^e_{xx} \frac{k_z^2}{2\omega\epsilon^2}
\\
-\sigma^m_{yy} + \frac{2k_z^1}{\omega\epsilon^1}
&
-\sigma^m_{yy} + \frac{2k_z^2}{\omega\epsilon^2}
\end{pmatrix}
\begin{pmatrix}
H_y^{(1,+)}
\\
H_y^{(2,-)}
\end{pmatrix}
\end{align}
```
"""

# ╔═╡ 860dc600-c86d-43a2-a88a-eb876eecbd6e
md"""
## Recovering Fresnel scattering

We can recover Fresnel scattering by setting the metasurface conductivities to zero, inverting the left hand side and reading off the scattered amplitudes.
To illustrate this, let's calculate this in this special case, which is just requiring the continuity of the tangential fields
```math
\begin{pmatrix}
I & -I
\\
\epsilon^2 k_z^1 & \epsilon^1 k_z^2
\end{pmatrix}
\begin{pmatrix}
H_y^{(1,-)}
\\
H_y^{(2,+)}
\end{pmatrix}
=
\begin{pmatrix}
-I & I
\\
\epsilon^2 k_z^1
& \epsilon^1 k_z^2
\end{pmatrix}
\begin{pmatrix}
H_y^{(1,+)}
\\
H_y^{(2,-)}
\end{pmatrix}
```
Let's invert the left hand side
```math
\begin{align}
\begin{pmatrix}
H_y^{(1,-)}
\\
H_y^{(2,+)}
\end{pmatrix}
&=
\begin{pmatrix}
(1 + \frac{\epsilon^2 k_z^1}{\epsilon^1 k_z^2})^{-1}
&
(\epsilon^1 k_z^2 + \epsilon^1 k_z^2)^{-1}
\\
-(1 + \frac{\epsilon^1 k_z^2}{\epsilon^2 k_z^1})^{-1}
&
(\epsilon^1 k_z^2 + \epsilon^2 k_z^1)^{-1}
\end{pmatrix}
\begin{pmatrix}
-I & I
\\
\epsilon^2 k_z^1 & \epsilon^1 k_z^2
\end{pmatrix}
\begin{pmatrix}
H_y^{(1,+)}
\\
H_y^{(2,-)}
\end{pmatrix}
\\ &=
(\epsilon^1 k_z^2 + \epsilon^2 k_z^1)^{-1}
\begin{pmatrix}
\epsilon^2 k_z^1 - \epsilon^1 k_z^2 & 2 \epsilon^1 k_z^2
\\
2 \epsilon^2 k_z^1 & \epsilon^1 k_z^2 - \epsilon^2 k_z^1
\end{pmatrix}
\begin{pmatrix}
H_y^{(1,+)}
\\
H_y^{(2,-)}
\end{pmatrix}
\end{align}
```
Note how the scattering depends only on the angle of incidence from the normal axis (the ratio of ``\omega`` and ``k_z``) and the impedances of the materials.
An analogous calculation can be done for the TE mode
"""

# ╔═╡ 271422af-db75-4373-be1c-104745f941bf
md"""
## Scattering matrix in 3d

The planewave form of Ampère's law we will use within dielectric medium ``i`` with permittivity ``\epsilon^i`` and permeability ``\mu^i`` is
```math
\bm{R} \bm{\Omega}^i \bm{R}^{-1} \bm{E}_\parallel^i = \omega \bm{R} \mu^i k_z^i \bm{H}_\parallel^i
```
where
```math
\bm{\Omega}^i = \bm{K} - \omega^2 \epsilon^i \mu^i \bm{I}
, \qquad
\bm{K} =
\begin{pmatrix}
k_x k_x & k_x k_y
\\
k_y k_x & k_y k_y
\end{pmatrix}
, \qquad
\bm{R} =
\begin{pmatrix}
0 & -1
\\
1 & 0
\end{pmatrix}
.
```
The GSTCs (with an ``\bm{R}`` and a 2 moved about) are
```math
\begin{align}
2 \bm{\rho}^e \bm{R} (\bm{H}_\parallel^2 - \bm{H}_\parallel^1)
&= \bm{E}_\parallel^2 + \bm{E}_\parallel^1
\\
\bm{E}_\parallel^2 - \bm{E}_\parallel^1
&= \bm{R} \bm{\sigma}^m (\bm{H}_\parallel^2 + \bm{H}_\parallel^1)/2
\end{align}
```
It is basic algebra to show that ``\bm{\Omega}^i`` commutes with ``\bm{\Omega}^j`` so multiplying the GSTCs by ``\bm{\Omega}^1 \bm{\Omega}^2 \bm{R}^{-1} = \bm{\Omega}^2 \bm{\Omega}^1 \bm{R}^{-1}`` gives
```math
\begin{align}
2 \bm{\Omega}^1 \bm{\Omega}^2 \bm{R}^{-1}
\bm{\rho}^e \bm{R} (\bm{H}_\parallel^2 - \bm{H}_\parallel^1)
&= \bm{\Omega}^1 \bm{R}^{-1} \bm{R} \bm{\Omega}^2 \bm{R}^{-1} \bm{E}_\parallel^2
+ \bm{\Omega}^2 \bm{R}^{-1} \bm{R} \bm{\Omega}^1 \bm{R}^{-1} \bm{E}_\parallel^1
\\
\bm{\Omega}^1 \bm{R}^{-1} \bm{R} \bm{\Omega}^2 \bm{R}^{-1} \bm{E}_\parallel^2
- \bm{\Omega}^2 \bm{R}^{-1} \bm{R} \bm{\Omega}^1 \bm{R}^{-1} \bm{E}_\parallel^1
&= \bm{\Omega}^1 \bm{\Omega}^2 \bm{\sigma}^m (\bm{H}_\parallel^2 + \bm{H}_\parallel^1)/2
\end{align}
```
Next we will do 3 operations in the same step:
1. Break each mode into forward and backward propagating components: ``\bm{H}_\parallel^i \to \bm{H}_\parallel^{(i,+)} + \bm{H}_\parallel^{(i,-)},\ \bm{E}_\parallel^i \to \bm{E}_\parallel^{(i,+)} + \bm{E}_\parallel^{(i,-)}``
2. Eliminate the electric field with Ampère's law: ``\bm{R} \bm{\Omega}^i \bm{R}^{-1} \bm{E}_\parallel^{(i,\pm)} = \pm \omega \bm{R} \mu^i k_z^i \bm{H}_\parallel^{(i,\pm)}``
3. Cancel any pairs of ``\bm{R}`` and ``\bm{R}^{-1}``
```math
\begin{align}
&2 \bm{\Omega}^1 \bm{\Omega}^2 \bm{R}^{-1} \bm{\rho}^e \bm{R} \left( \bm{H}_\parallel^{(2,+)} + \bm{H}_\parallel^{(2,-)} - \bm{H}_\parallel^{(1,+)} - \bm{H}_\parallel^{(1,-)} \right)
\\ &\qquad = \omega \left(
\bm{\Omega}^1 \mu^2 k_z^2 \left( \bm{H}_\parallel^{(2,+)} - \bm{H}_\parallel^{(2,-)} \right)
+ \bm{\Omega}^2 \mu^1 k_z^1 \left( \bm{H}_\parallel^{(1,+)} - \bm{H}_\parallel^{(1,-)} \right) \right)
\\
&\omega \left(
\bm{\Omega}^1 \mu^2 k_z^2 \left( \bm{H}_\parallel^{(2,+)} - \bm{H}_\parallel^{(2,-)} \right)
- \bm{\Omega}^2 \mu^1 k_z^1 \left( \bm{H}_\parallel^{(1,+)} - \bm{H}_\parallel^{(1,-)} \right) \right)
\\ &\qquad = \bm{\Omega}^1 \bm{\Omega}^2 \bm{\sigma}^m \left( \bm{H}_\parallel^{(2,+)} + \bm{H}_\parallel^{(2,-)} + \bm{H}_\parallel^{(1,+)} + \bm{H}_\parallel^{(1,-)} \right)/2
\end{align}
```
Next if we move the scattered components, (1,-) and (2,+), to the left and the incident components to the right, (1,+) and (2,-), we can express the scattering matrix as
```math
\begin{align}
&\begin{pmatrix}
-\bm{\Omega}^1 \bm{\Omega}^2 \bm{R}^{-1} \bm{\rho}^e \bm{R} + \frac{1}{2} \omega \bm{\Omega}^2 \mu^1 k_z^1
&
\bm{\Omega}^1 \bm{\Omega}^2 \bm{R}^{-1} \bm{\rho}^e \bm{R} - \frac{1}{2} \omega \bm{\Omega}^1 \mu^2 k_z^2
\\
-\bm{\Omega}^1 \bm{\Omega}^2 \bm{\sigma}^m + 2 \omega \bm{\Omega}^2 \mu^1 k_z^1
&
-\bm{\Omega}^1 \bm{\Omega}^2 \bm{\sigma}^m + 2 \omega \bm{\Omega}^1 \mu^2 k_z^2
\end{pmatrix}
\begin{pmatrix}
\bm{H}_\parallel^{(1,-)}
\\
\bm{H}_\parallel^{(2,+)}
\end{pmatrix}
\\ =
&\begin{pmatrix}
\bm{\Omega}^1 \bm{\Omega}^2 \bm{R}^{-1} \bm{\rho}^e \bm{R} + \frac{1}{2} \omega \bm{\Omega}^2 \mu^1 k_z^1
&
-\bm{\Omega}^1 \bm{\Omega}^2 \bm{R}^{-1} \bm{\rho}^e \bm{R} - \frac{1}{2} \omega \bm{\Omega}^1 \mu^2 k_z^2
\\
\bm{\Omega}^1 \bm{\Omega}^2 \bm{\sigma}^m + 2 \omega \bm{\Omega}^2 \mu^1 k_z^1
&
\bm{\Omega}^1 \bm{\Omega}^2 \bm{\sigma}^m + 2 \omega \bm{\Omega}^1 \mu^2 k_z^2
\end{pmatrix}
\begin{pmatrix}
\bm{H}_\parallel^{(1,+)}
\\
\bm{H}_\parallel^{(2,-)}
\end{pmatrix}
\end{align}
```
Note that if the two media are the same (i.e. ``\epsilon^1=\epsilon^2, \mu^1=\mu^2``) then ``\bm{\Omega}^1 = \bm{\Omega}^2`` and the scattering matrix has an overall factor of ``\Omega`` that can be removed, which recovers the result for the scattering matrix in a homogenous medium.
In fact we can manipulate this expression further using the specific form of the inverse
```math
(\bm{\Omega}^i)^{-1} = \bm{R}^{-1} \bm{\Omega}^i \bm{R} \frac{1}{\omega^2 \epsilon^i \mu^i (k_z^i)^2}
```
due to the fact ``\bm{\Omega}^i`` is a square block matrix whose blocks are diagonal matrices.
"""

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

Kuester published a follow-up article (DOI: 10.1109/APS.2010.5562250)
with a first order correction to the previous one which should generalize to the
case where the metafilm has different media on either side.
The form of the GSTC is the same as before, except the impedance matrices are modified functions of the actual material parameters
"

# ╔═╡ Cell order:
# ╟─83d9578c-3960-486c-ad3d-b4449f7bae6f
# ╠═c24044bc-2fc3-4e34-b99c-7dc4ec6a7515
# ╠═42da96ce-8122-4aa0-99d6-635f68db6e66
# ╠═88e46392-1fbc-4079-8207-dd5b3ac01d7d
# ╠═6cbc5f34-fe14-45e8-9fee-9787d76ad0fe
# ╠═fb8df44c-a354-42de-9e99-df64dd4c5727
# ╠═f3d36c35-c09b-483f-bb41-046b6ecd0773
# ╠═e899dbc0-03bb-4a13-9c09-bfad33e5c1d3
# ╠═43d29086-510d-4870-94b7-2ff5ff34324b
# ╠═5377ae36-32fc-41c5-bd03-c737422b17f4
# ╠═65c2e966-8cb0-4d76-96db-18f9045695a3
# ╠═16316725-c608-49a5-b90c-8caff973683a
# ╠═bf1213cd-2ec1-4d1e-8ebc-5be5af26621c
# ╟─ba0c4060-f9d4-43b6-bc31-b48dcd81dd91
# ╟─687627f9-6938-4d7a-9395-9f1783b8a449
# ╟─860dc600-c86d-43a2-a88a-eb876eecbd6e
# ╟─271422af-db75-4373-be1c-104745f941bf
# ╟─325b8359-ebac-434b-bba1-31aba05c2c39
