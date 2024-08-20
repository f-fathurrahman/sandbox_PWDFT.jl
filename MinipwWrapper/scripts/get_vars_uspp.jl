using Serialization: serialize

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )

function get_vars_uspp()

    nhm = unsafe_load(cglobal((:__uspp_param_MOD_nhm, LIBMINIPW), Int32)) |> Int64
    Natoms = unsafe_load(cglobal((:__ions_base_MOD_nat, LIBMINIPW), Int32)) |> Int64

    # qq_at(nhm,nhm,nat)
    ptr = cglobal((:__uspp_MOD_qq_at, LIBMINIPW), Ptr{Float64})
    Ndim1 = nhm
    Ndim2 = nhm
    Ndim3 = Natoms
    tmp = zeros(Float64, Ndim1*Ndim2*Ndim3)
    ip = 1
    for l in 1:Ndim3, j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    ddd_paw = reshape(tmp, Ndim1, Ndim2, Ndim3)
    serialize("qq_at.dat", ddd_paw)

    # Hardcoded parameters in uspp module
    lmaxx  = 3
    lqmax = 2*lmaxx + 1
    nlx = (lmaxx + 1)^2
    mx = 2*lqmax-1

    # uspp_ap
    # ap is an automatic array, can simply use unsafe_wrap
    Ndim1 = lqmax^2
    Ndim2 = nlx
    Ndim3 = nlx
    ptr = cglobal( (:__uspp_MOD_ap, LIBMINIPW), Float64 )
    uspp_ap = unsafe_wrap(Array{Float64,3}, ptr, (Ndim1,Ndim2,Ndim3))
    serialize("uspp_ap.dat", uspp_ap)

    # lpx(nlx,nlx) ! maximum combined angular momentum LM
    Ndim1 = nlx
    Ndim2 = nlx
    ptr = cglobal( (:__uspp_MOD_lpx, LIBMINIPW), Int32 )
    uspp_lpx = unsafe_wrap(Array{Int32,2}, ptr, (Ndim1,Ndim2))
    serialize("uspp_lpx.dat", uspp_lpx)


    # lpl(nlx,nlx,mx) ! list of combined angular momenta  LM
    Ndim1 = nlx
    Ndim2 = nlx
    Ndim3 = mx
    ptr = cglobal( (:__uspp_MOD_lpl, LIBMINIPW), Int32 )
    uspp_lpl = unsafe_wrap(Array{Int32,3}, ptr, (Ndim1,Ndim2,Ndim3))
    serialize("uspp_lpl.dat", uspp_lpl)

end

#=
ap = deserialize("uspp_ap.dat");
lpx = convert(Array{Int64,2}, deserialize("uspp_lpx.dat"));
lpl = convert(Array{Int64,3}, deserialize("uspp_lpl.dat"));
=#


get_vars_uspp()

