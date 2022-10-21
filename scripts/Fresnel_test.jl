
ω = 8.0
kx = 1.0
# vacuum
ϵ₁ = 1.0
μ₁ = 1.0
# water
ϵ₂ = 1.77
μ₂ = 1.0

S = uFresnelTM(ω, kx, ϵ₁, μ₁, ϵ₂, μ₂)
J = [
    0 1
    1 0
]

sheet = TransparentSheet()

n = 5
L = 1.0
# 2D
dims = ((n, L), )
nrep = 2
# 3D
# dims = ((n, L), (n,L))
# nrep = 4
pw = PlaneWaves(ω, dims)
M1 = UniformMedium{ϵ₁,μ₁}()
M2 = UniformMedium{ϵ₂,μ₂}()
DS = smatrix(HField(), pw, sheet, M1, M2)

BDI = map(i -> CartesianIndex(blockify.(i.I, (n,n))), CartesianIndices(DS))
BDS = similar(DS)
for (i, e) in zip(BDI, DS)
    BDS[i] = e
end

k₁ = ω*sqrt(ϵ₁*μ₁)
k₂ = ω*sqrt(ϵ₂*μ₂)

test = [
    (ϵ₂*k₁ - ϵ₁*k₂)   2ϵ₂*k₁
    2ϵ₁*k₂   -(ϵ₂*k₁ - ϵ₁*k₂)
] ./ (ϵ₂*k₁ + ϵ₁*k₂)