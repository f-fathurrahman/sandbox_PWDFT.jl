# Test naive FFT

using LinearAlgebra: norm
using PWDFT


include("GVectorsFull.jl")

function func_R_analytic(rl::Float64; rloc=1.5)
    return exp(-0.5*(rl/rloc)^2)
end

function test_eval_R( LatVecs, func_R, r, offset_neighbors; verbose=false )
    # real space evalution
    rl = norm(r) # should be distance
    fR = func_R(rl) # This should be enough if the function is localized in unit cell
    #
    v1 = LatVecs[:,1]
    v2 = LatVecs[:,2]
    v3 = LatVecs[:,3]
    nnx = offset_neighbors[1]; @assert nnx >= 0
    nny = offset_neighbors[2]; @assert nny >= 0
    nnz = offset_neighbors[3]; @assert nnz >= 0
    num_neighbors = 0
    for i in -nnx:nnx, j in -nny:nny, k in -nnz:nnz
        if i == 0 && j == 0 && k == 0
            continue
        end
        num_neighbors += 1
        rn = r + i*v1 + j*v2 + k*v3
        rl = norm(rn)
        fR += func_R(rl)
    end
    verbose && @info "num_neighbors = $(num_neighbors)"
    return fR
end

# with G-vectors full, given G-vector index ig
function naive_R_to_G_with_gvec_full(pw, fR, ig)
    
    Nx, Ny, Nz = pw.Ns
    fRr = reshape(fR, Nx, Ny, Nz)
    fG = zeros(ComplexF64, Nx, Ny, Nz)
    #
    gvec_full = GVectorsFull(pw.Ns, pw.RecVecs)
    j, k, l = gvec_full.idx_g2miller[ig]
    @info "FFT grid index ip = $(gvec_full.idx_g2r[ig])"
    @info "Original $j, $k, $l"
    #
    if j < 0
        j = j + Nx
    end
    if k < 0
        k = k + Ny
    end
    if l < 0
        l = l + Nz
    end
    @info "After adding offset: $j $k $l"
    #
    res = 0.0 + im*0.0
    for w in 0:(Nz-1), v in 0:(Ny-1), u in 0:(Nx-1)
        fex = exp(-im*2*pi/Nx*j*u)
        fey = exp(-im*2*pi/Ny*k*v)
        fez = exp(-im*2*pi/Nz*l*w)
        res += fRr[u+1,v+1,w+1] * fex * fey * fez
        # need to offset the indices u,v,w with 1
    end
    println("res = ", res)

end

function naive_R_to_G(pw, fR)
    #
    Nx, Ny, Nz = pw.Ns
    fRr = reshape(fR, Nx, Ny, Nz) # input
    fG = zeros(ComplexF64, Nx, Ny, Nz) # output
    #
    for l in 0:(Nz-1), k in 0:(Ny-1), j in 0:(Nx-1)
        res = 0.0 + im*0.0
        for w in 0:(Nz-1), v in 0:(Ny-1), u in 0:(Nx-1)
            fex = exp(-im*2*pi/Nx*j*u)
            fey = exp(-im*2*pi/Ny*k*v)
            fez = exp(-im*2*pi/Nz*l*w)
            res += fRr[u+1,v+1,w+1] * fex * fey * fez
            # need to offset the indices u,v,w with 1
        end
        fG[j+1,k+1,l+1] = res # offset indices j,k,l with 1
    end
    return reshape(fG, prod(pw.Ns))
end



L = 10.0 # bohr
ecutwfc = 5.0
LatVecs = gen_lattice_fcc(L)
rloc = 1.4
func_R(r) = func_R_analytic(r, rloc=rloc)
pw = PWGrid(ecutwfc, LatVecs)
Npoints = prod(pw.Ns)
fR_grid = zeros(Float64, Npoints)
r_grid = PWDFT.init_grid_R(pw.Ns, LatVecs)

offset_neighbors = [1, 1, 1]
for ip in 1:Npoints
    fR_grid[ip] = test_eval_R(LatVecs, func_R, r_grid[:,ip], offset_neighbors)
end

ctmp = zeros(ComplexF64, Npoints)
ctmp[:] .= fR_grid[:]
R_to_G!(pw, ctmp)
# scale with 1/Npoints
#ctmp[:] *= (1/Npoints)