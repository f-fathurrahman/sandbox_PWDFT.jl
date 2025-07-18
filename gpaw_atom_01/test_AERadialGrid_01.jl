using Printf

include("AERadialGrid.jl")

maxnodes = 2
gpernode = 150 # default?
β = 0.4 # parameter for radial grid
N = gpernode*(maxnodes + 1)

a = β/N
b = 1/N
ae_grid = AERadialGrid(a, b, npts=N)
for i in 1:N
    @printf("%4d %18.10f %18.10f %18.10f\n", i,
        ae_grid.r[i], ae_grid.dr[i], ae_grid.d2gdr2[i])
end
