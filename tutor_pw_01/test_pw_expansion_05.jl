# Like test_pw_expansion_03, but with structure factor

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

function test_eval_G(LatVecs, ecutwfc, func_G, r, r0)
    pw = PWGrid(ecutwfc, LatVecs)
    Ng = pw.gvec.Ng
    #@info "Ng = $Ng"
    G2 = pw.gvec.G2
    G = pw.gvec.G
    CellVolume = pw.CellVolume
    fG = zeros(Float64, Ng)
    # Evaluate the coefficients
    for ig in 1:Ng
        Gl = sqrt(G2[ig])
        fG[ig] = func_G(Gl)/CellVolume # we include 1/CellVolume factor here
    end
    Sf = zeros(ComplexF64, Ng)
    for ig in 1:Ng
        Sf[ig] = exp(-im*dot(G[:,ig], r0))
    end
    #
    s = 0.0 + im*0.0
    for ig in 1:Ng
        Gr = dot(G[:,ig], r)
        s += fG[ig] * exp(im*Gr) * Sf[ig]
    end
    return s
end


function test_eval_R( LatVecs, func_R, r, r0, offset_neighbors )
    # real space evalution
    rl = norm(r - r0) # should be distance
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
        rl = norm(rn - r0) # don't forget to subtract with center
        fR += func_R(rl)
    end
    #@info "num_neighbors = $(num_neighbors)"
    return fR
end


function main_recip_space(r, r0, rloc)
    println("\nTest evaluating function in reciprocal space, vary ecutwfc")
    L = 10.0 # bohr
    LatVecs = gen_lattice_sc(L)
    println("rloc = $rloc, smaller rloc will need larger cutoff")
    func_G(G) = func_G_analytic(G, rloc=rloc) # set rloc
    #
    results = ComplexF64[]
    ecutwfc_list = 1.0:2.0:10.0
    for ecutwfc in ecutwfc_list
        res_G = test_eval_G(LatVecs, ecutwfc, func_G, r, r0)# |> real
        append!(results, res_G)
    end
    ref = results[end]
    for i in 1:length(ecutwfc_list)
        # Take the last for the converged value
        println("ecutwfc = $(ecutwfc_list[i]) res = $(results[i]) abs diff = $(abs(results[i]-ref))")
    end
    # More localized functions (smaller rloc) will need larger ecut to converge
    #
    return results[end]
end


function main_real_space(r, r0, rloc)
    #
    println("\nTest evaluating function in real space, vary num_neighbors")
    println("rloc = $rloc, larger rloc (more delocalized) will need larger num_neighbors")
    func_R(r) = func_R_analytic(r, rloc=rloc) # set rloc
    results = Float64[]
    nn_list = 0:1:10
    for nn in nn_list
        res_R = test_eval_R(LatVecs, func_R, r, r0, [nn,nn,nn])
        append!(results, res_R)
    end
    ref = results[end]
    for i in 1:length(nn_list)
        # Take the last nn for the converged value
        println("nn = $(nn_list[i]) , res = $(results[i]) abs diff = $(abs(results[i]-ref))")
    end
    # More delocalized functions will need larger number of neighbors to converge
    println("Converged result = ", results[end])
    #
    return results[end]
end


function test_main()
    # Must use the same rloc, r, and r0
    rloc = 1.0
    r = [5.4, 5.4, 5.4] # point of evaluation
    r0 = [5.0, 5.1, 6.0]
    fG = main_recip_space(r, r0, rloc)
    fR = main_real_space(r, r0, rloc)
    println("fG = $fG")
    println("fR = $fR")
end

test_main()

