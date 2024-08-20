using Serialization: serialize

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )

function get_ddd_paw()

    nhm = unsafe_load(cglobal((:__uspp_param_MOD_nhm, LIBMINIPW), Int32)) |> Int64
    Natoms = unsafe_load(cglobal((:__ions_base_MOD_nat, LIBMINIPW), Int32)) |> Int64
    Nspin = unsafe_load(cglobal((:__lsda_mod_MOD_nspin, LIBMINIPW), Int32)) |> Int64

    # ddd_paw(nhm*(nhm+1)/2,nat,nspin) )
    ptr = cglobal((:__paw_variables_MOD_ddd_paw, LIBMINIPW), Ptr{Float64})
    Ndim1 = round(Int64, nhm*(nhm+1)/2)
    Ndim2 = Natoms
    Ndim3 = Nspin
    tmp = zeros(Float64, Ndim1*Ndim2*Ndim3)
    ip = 1
    for l in 1:Ndim3, j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    ddd_paw = reshape(tmp, Ndim1, Ndim2, Ndim3)
    serialize("ddd_paw.dat", ddd_paw)

end

get_ddd_paw()