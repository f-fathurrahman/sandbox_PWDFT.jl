using Serialization: serialize

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )

function pwx_uspp_param()

    nsp = unsafe_load(cglobal((:__ions_base_MOD_nsp, LIBMINIPW), Int32)) |> Int64
    Nspecies = nsp
    
    ptr = cglobal((:__uspp_param_MOD_nh, LIBMINIPW), Int32)
    nh = zeros(Int64,nsp)
    for isp in 1:nsp
        nh[isp] = unsafe_load(ptr,isp) |> Int64
    end

    nhm = unsafe_load(cglobal((:__uspp_param_MOD_nhm, LIBMINIPW), Int32)) |> Int64
    # nhm is max(nh(1:ntyp))
    println("nsp = ", nsp)
    
    println("nh = ", nh)
    println("nhm = ", nhm)
    # reference: init_run.f90 subroutine pre_init()


    # indv, nhtol, nhtolm: allocated in allocate_nlpot
    # set in init_us_1
    # defined in uspp module

    ptr = cglobal((:__uspp_MOD_indv, LIBMINIPW), Ptr{Int32})
    Ndim1 = nhm
    Ndim2 = Nspecies
    tmp = zeros(Int64, Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip) |> Int64
        ip = ip + 1
    end
    indv = reshape(tmp, nhm, Nspecies)
    println("indv = ")
    for isp in 1:Nspecies
        println(indv[1:nh[isp],isp])
    end

    ptr = cglobal((:__uspp_MOD_nhtol, LIBMINIPW), Ptr{Int32})
    Ndim1 = nhm
    Ndim2 = Nspecies
    tmp = zeros(Int64, Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip) |> Int64
        ip = ip + 1
    end
    nhtol = reshape(tmp, nhm, Nspecies)
    println("nhtol = ")
    for isp in 1:Nspecies
        println(nhtol[1:nh[isp],isp])
    end


    ptr = cglobal((:__uspp_MOD_nhtolm, LIBMINIPW), Ptr{Int32})
    Ndim1 = nhm
    Ndim2 = Nspecies
    tmp = zeros(Int64, Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip) |> Int64
        ip = ip + 1
    end
    nhtolm = reshape(tmp, nhm, Nspecies)
    println("nhtolm = ")
    for isp in 1:Nspecies
        println(nhtolm[1:nh[isp],isp])
    end

    #indv = zeros(Int64, nhm, Nspecies)
    #nhtol = zeros(Int64, nhm, Nspecies)
    #nhtolm = zeros(Int64, nhm, Nspecies)

end

pwx_uspp_param()