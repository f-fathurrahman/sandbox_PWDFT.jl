function elk_write_atomic_species_scalars()
    nspecies = unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
end

function elk_write_nrsp()
    nspecies = unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
    nrsp = zeros(Int64,nspecies)
    ptr = cglobal((:__m_atomic_species_MOD_nrsp, LIBLAPW), Int32)
    for i in 1:nspecies
        nrsp[i] = unsafe_load(ptr,i)
    end
    serialize(joinpath(ELK_DATA_DIR,"nrsp.dat"), nrsp)
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
    serialize(joinpath(ELK_DATA_DIR,"rsp.dat"), rsp)
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
    serialize(joinpath(ELK_DATA_DIR,"rhosp.dat"), rhosp)
    return
end