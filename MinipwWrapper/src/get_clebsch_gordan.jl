using Serialization: serialize

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )

function pwx_get_clebsch_gordan()
    # Hardcoded parameters
    lmaxx  = 3
    lqmax = 2*lmaxx + 1
    nlx = (lmaxx + 1)^2

    # ap is automatic array
    Ndim1 = lqmax^2
    Ndim2 = nlx
    Ndim3 = nlx
    ptr = cglobal( (:__uspp_MOD_ap, LIBMINIPW), Float64 )
    uspp_ap = unsafe_wrap(Array{Float64,3}, ptr, (Ndim1,Ndim2,Ndim3))
    serialize("uspp_ap.dat", uspp_ap)
end

pwx_get_clebsch_gordan()
