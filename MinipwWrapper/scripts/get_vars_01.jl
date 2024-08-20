using Serialization: serialize

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )


function pwx_get_evc()
    npwx = unsafe_load(cglobal((:__wvfct_MOD_npwx, LIBMINIPW), Int32)) |> Int64
    nbnd = unsafe_load(cglobal((:__wvfct_MOD_nbnd, LIBMINIPW), Int32)) |> Int64
    npol = unsafe_load(cglobal((:__noncollin_module_MOD_npol, LIBMINIPW), Int32)) |> Int64
    

    nks = unsafe_load(cglobal((:__klist_MOD_nks, LIBMINIPW), Int32)) |> Int64
    println("nks = ", nks)
    # Read ngk
    ptr = cglobal( (:__klist_MOD_ngk, LIBMINIPW), Ptr{Int32} )
    ngk = zeros(Int64, nks)
    for i in 1:nks
        ngk[i] = unsafe_load(unsafe_load(ptr,1),i) |> Int64
    end
    println("ngk = ", ngk)


    ptr = cglobal( (:__wavefunctions_MOD_evc, LIBMINIPW), Ptr{ComplexF64} )
    println("npwx = ", npwx)
    println("nbnd = ", nbnd)
    println("npol = ", npol)
    Ndim1 = npwx*npol
    Ndim2 = nbnd
    tmp = zeros(ComplexF64,npwx*npol*nbnd)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    evc = reshape(tmp, Ndim1, Ndim2)
    # evc(npwx*npol,nbnd)
    serialize("evc.dat", evc)
    println("Pass here")
end

pwx_get_evc()
