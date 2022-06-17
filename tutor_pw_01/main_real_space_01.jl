using LinearAlgebra: diagm
using PWDFT

#LatVecs = gen_lattice_fcc(16.0)
LatVecs = gen_lattice_tetragonal_P(10.0, 12.0)
display(LatVecs); println()

h = LatVecs
a1 = LatVecs[:,1]
a2 = LatVecs[:,2]
a3 = LatVecs[:,3]

display(a1); println()
display(a2); println()
display(a3); println()

Nx = 2
Ny = 3
Nz = 5
NN = diagm(0 => [1/Nx, 1/Ny, 1/Nz])

qx = range(0, Nx, step=1)
qy = range(0, Ny, step=1)
qz = range(0, Nz, step=1)

Npoints = Nx * Ny * Nz
R = zeros(Float64, 3, Npoints)
ip = 0
q = zeros(Float64, 3)
for l in 1:Nz, k in 1:Ny, j in 1:Nx
    q[1] = qx[j]
    q[2] = qy[k]
    q[3] = qz[l]
    ip = ip + 1
    R[:,ip] = h*NN*q
end

display(R'); println()

