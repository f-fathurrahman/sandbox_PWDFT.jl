using Printf
using LinearAlgebra

using PWDFT

function eval_func_G(Gl::Float64; rloc=1.0)
    pre2 = sqrt(8*π^3)*rloc^3
    expGr2 = exp(-0.5*(Gl*rloc)^2)
    return pre2*expGr2
end

function eval_func_R(rl::Float64; rloc=1.0)
    return exp(-0.5*(rl/rloc)^2)
end

L = 10.0 # bohr
ecutwfc = 5.0 # hartree
# Using higher cutoff will give better results

LatVecs = gen_lattice_sc(L)
pw = PWGrid(ecutwfc, LatVecs)

println(pw)

Ng = pw.gvec.Ng
G2 = pw.gvec.G2
G = pw.gvec.G
CellVolume = pw.CellVolume

fG = zeros(Float64, Ng)
for ig in 1:Ng
    Gl = sqrt(G2[ig])
    fG[ig] = eval_func_G(Gl)/CellVolume
end


Npoints = prod(pw.Ns)
Rgrid = PWDFT.init_grid_R(pw.Ns, LatVecs)

println("Evaluating plane wave expansion")

ip = 10
r = Rgrid[:,ip]
s = 0.0 + im*0.0
for ig in 1:Ng
    Gr = dot(G[:,ig], r)
    s += fG[ig] * exp(im*Gr)
end

fR = real(s)
println("PW expansion = ", fR)

rl = norm(r)
fR_a = eval_func_R(rl)
println("analytic = ", fR_a)
println("Difference = ", abs(fR - fR_a))

v1 = LatVecs[:,1]
v2 = LatVecs[:,2]
v3 = LatVecs[:,3]

# Add nearest neighbor contributions
# Result will depend on rloc and CellSize
NumNeighbors = 0
for i in -1:1, j in -1:1, k in -1:1
    if i == 0 && j == 0 && k == 0
        continue
    end
    NumNeighbors += 1
    rn = r + i*v1 + j*v2 + k*v3
    rl = norm(rn)
    fR_a += eval_func_R(rl)
end
println("NumNeighbors = ", NumNeighbors) # should be 27 - 1 = 26

println("analytic (adding nearest neighbors contrib) = ", fR_a)
println("Difference = ", abs(fR - fR_a))
