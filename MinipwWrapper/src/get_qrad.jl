using Serialization: serialize

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )

function pwx_get_qrad()

    nsp = unsafe_load(cglobal((:__ions_base_MOD_nsp, LIBMINIPW), Int32)) |> Int64
    nbetam = unsafe_load(cglobal((:__uspp_param_MOD_nbetam, LIBMINIPW), Int32)) |> Int64
    lmaxq = unsafe_load(cglobal((:__uspp_param_MOD_lmaxq, LIBMINIPW), Int32)) |> Int64
    nqxq = unsafe_load(cglobal((:__us_MOD_nqxq, LIBMINIPW), Int32)) |> Int64
    
    println("nsp = ", nsp)
    println("nbetam = ", nbetam)
    println("lmaxq = ", lmaxq)
    println("nqxq = ", nqxq)

    ptr = cglobal( (:__us_MOD_qrad, LIBMINIPW), Ptr{Float64} )
    Ndim1 = nqxq
    Ndim2 = Int64(nbetam*(nbetam+1)/2)
    Ndim3 = lmaxq
    Ndim4 = nsp
    tmp = zeros(Float64, Ndim1*Ndim2*Ndim3*Ndim4)
    ip = 1
    for l in 1:Ndim4, k in 1:Ndim3, j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    qrad = reshape(tmp, Ndim1, Ndim2, Ndim3, Ndim4)
    serialize("qrad.dat", qrad)
    #qrad(nqxq, nbetam*(nbetam+1)/2, lmaxq, nsp)
    return
end

pwx_get_qrad()