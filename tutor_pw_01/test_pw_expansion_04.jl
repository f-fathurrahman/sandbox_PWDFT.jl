using Printf
using LinearAlgebra: dot, norm
using PWDFT

function func_G_analytic(Gl::Float64; rloc=1.5)
    pre2 = sqrt(8*π^3)*rloc^3
    expGr2 = exp(-0.5*(Gl*rloc)^2)
    return pre2*expGr2
end

function func_R_analytic(rl::Float64; rloc=1.5)
    return exp(-0.5*(rl/rloc)^2)
end

function test_eval_G(LatVecs, ecutwfc, func_G, r; verbose=false)
    pw = PWGrid(ecutwfc, LatVecs)
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    G = pw.gvec.G
    CellVolume = pw.CellVolume
    fG = zeros(Float64, Ng)
    # Evaluate the coefficients
    for ig in 1:Ng
        Gl = sqrt(G2[ig])
        fG[ig] = func_G(Gl)/CellVolume # we include 1/CellVolume factor here
    end
    s = 0.0 + im*0.0
    for ig in 1:Ng
        Gr = dot(G[:,ig], r)
        s += fG[ig] * exp(im*Gr)
    end
    verbose && @info "Ng = $Ng"
    return s
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



#function main()

    L = 10.0 # bohr
    ecutwfc = 5.0
    LatVecs = gen_lattice_sc(L)
    rloc = 1.4
    func_G(G) = func_G_analytic(G, rloc=rloc)
    func_R(r) = func_R_analytic(r, rloc=rloc)
    
    pw = PWGrid(ecutwfc, LatVecs)
    Npoints = prod(pw.Ns)
    fR_grid = zeros(Float64, Npoints)
    r_grid = PWDFT.init_grid_R(pw.Ns, LatVecs)
    
    #=
    # Testing nn convergence
    r = r_grid[:,10]
    @info "r = $r"
    results = Float64[]
    for nn in [0, 1]
        res_R = test_eval_R(LatVecs, func_R, r, [nn,nn,nn])
        append!(results, res_R)
    end
    println("Differences = ", abs.(results[end] .- results))
    # More delocalized functions will need larger number of neighbors to converge
    =#

    offset_neighbors = [1, 1, 1]
    for ip in 1:Npoints
        fR_grid[ip] = test_eval_R(LatVecs, func_R, r_grid[:,ip], offset_neighbors)
    end
    dVol = pw.CellVolume/Npoints # evaluate its integral within unit cell
    println("integ fR_grid = ", sum(fR_grid)*dVol)
    # println("Analytic integ = ", sqrt(2*pi)*rloc) # not applicable because of its neighbors

    # We now evaluate the expansion coefs using FFT
    ctmp = zeros(ComplexF64, Npoints)
    ctmp[:] .= fR_grid[:]
    R_to_G!(pw, ctmp)
    # scale with 1/Npoints
    ctmp[:] *= (1/Npoints)
    #
    fG = zeros(ComplexF64, pw.gvec.Ng)
    for ig in 1:Ng
        ip = pw.gvec.idx_g2r[ig]
        fG[ig] = ctmp[ip]
    end
    # Compare with analytic
    fG_analytic = zeros(ComplexF64, pw.gvec.Ng)
    for ig in 1:Ng
        Gl = sqrt(pw.gvec.G2[ig])
        fG_analytic[ig] = func_G(Gl)/pw.CellVolume # need to rescale with 1/CellVolume
    end

#end
#main()