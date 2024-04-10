function get_nrspmax()
    nrspmax = unsafe_load(cglobal((:__m_atomic_species_MOD_nrspmax, LIBLAPW), Int32)) |> Int64
    return nrspmax
end

# nrsp is an automatic array (with preallocated size of maxspecies)
function get_nrsp()
    symbol = :__m_atomic_species_MOD_nrsp
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end


# Radial grid for each species
function get_rsp()
    symbol = :__m_atomic_species_MOD_rsp
    nrspmax  = get_nrspmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (nrspmax,nspecies))
end


# Atomic rho
function get_rhosp()
    symbol = :__m_atomic_species_MOD_rhosp
    nrspmax  = get_nrspmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (nrspmax,nspecies))
end