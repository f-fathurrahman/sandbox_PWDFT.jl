using OffsetArrays
using Serialization: serialize
const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/src_LAPW/liblapwdft.so"

function elk_init_run()
    ccall( (:read_input_, LIBLAPW), Cvoid, () )
    ccall( (:init0_, LIBLAPW), Cvoid, () )
    ccall( (:init1_, LIBLAPW), Cvoid, () )
    ccall( (:info_gvectors_, LIBLAPW), Cvoid, () )
    ccall( (:info_muffin_tins_, LIBLAPW), Cvoid, () )
    ccall( (:writesym_, LIBLAPW), Cvoid, () )
    isdir("TEMP_datadir") || mkdir("TEMP_datadir")
    return
end


function elk_write_nrsp()
    nspecies = unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
    nrsp = zeros(Int64,nspecies)
    ptr = cglobal((:__m_atomic_species_MOD_nrsp, LIBLAPW), Int32)
    for i in 1:nspecies
        nrsp[i] = unsafe_load(ptr,i)
    end
    serialize("TEMP_datadir/nrsp.dat", nrsp)
    return
end

# Radial grid for each species
function elk_write_rsp()
    #
    nrspmax  = unsafe_load(cglobal((:__m_atomic_species_MOD_nrspmax, LIBLAPW), Int32)) |> Int64
    nspecies = unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
    #
    ptr_rsp = cglobal( (:__m_atomic_species_MOD_rsp, LIBLAPW), Ptr{Float64} )
    rsp = zeros(Float64, nrspmax*nspecies)
    ip = 1
    for j in 1:nspecies, i in 1:nrspmax
        rsp[ip] = unsafe_load(unsafe_load(ptr_rsp,1),ip)
        ip = ip + 1
    end
    rsp = reshape(rsp, (nrspmax,nspecies))
    serialize("TEMP_datadir/rsp.dat", rsp)
    return
end


# Atomic rho
function elk_write_rhosp()
    #
    nrspmax  = unsafe_load(cglobal((:__m_atomic_species_MOD_nrspmax, LIBLAPW), Int32)) |> Int64
    nspecies = unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
    #
    ptr_rhosp = cglobal( (:__m_atomic_species_MOD_rhosp, LIBLAPW), Ptr{Float64} )
    rhosp = zeros(Float64, nrspmax*nspecies)
    ip = 1
    for j in 1:nspecies, i in 1:nrspmax
        rhosp[ip] = unsafe_load(unsafe_load(ptr_rhosp,1),ip)
        ip = ip + 1
    end
    rhosp = reshape(rhosp, (nrspmax,nspecies))
    serialize("TEMP_datadir/rhosp.dat", rhosp)
    return
end


function elk_write_rlmt()
    #
    lmaxo = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxo, LIBLAPW), Int32 )) |> Int64
    nrmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_nrmtmax, LIBLAPW), Int32 )) |> Int64
    #
    Ndim1 = nrmtmax
    Ndim2 = 2*(lmaxo+2)
    Ndim3 = nspecies
    rlmt = zeros(Float64, Ndim1*Ndim2*Ndim3)
    ptr_rlmt = cglobal( (:__m_muffin_tins_MOD_rlmt, LIBLAPW), Ptr{Float64} )
    ip = 1
    for k in 1:Ndim3, j in 1:Ndim2, i in 1:Ndim1
        rlmt[ip] = unsafe_load(unsafe_load(ptr_rlmt,1),ip)
        ip = ip + 1
    end
    rlmt = reshape(rlmt, (Ndim1, Ndim2, Ndim3))
    rlmt = OffsetArray(rlmt, 1:nrmtmax, -lmaxo-1:lmaxo+2, 1:nspecies)
    serialize("TEMP_datadir/rlmt.dat")
end