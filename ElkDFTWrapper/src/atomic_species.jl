# nstsp is an automatic array (with preallocated size of maxspecies)
# only load nspecies data
function get_nstsp()
    symbol = :__m_atomic_species_MOD_nstsp
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end

function get_nstspmax()
    nstspmax = unsafe_load(cglobal((:__m_atomic_species_MOD_nstspmax, LIBLAPW), Int32)) |> Int64
    return nstspmax
end

function get_nrspmax()
    nrspmax = unsafe_load(cglobal((:__m_atomic_species_MOD_nrspmax, LIBLAPW), Int32)) |> Int64
    return nrspmax
end

# nrsp is an automatic array (with preallocated size of maxspecies)
# XXX: We only load 1:nspecies data
function get_nrsp()
    symbol = :__m_atomic_species_MOD_nrsp
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end

# effective infinity for species
function get_rmaxsp()
    symbol = :__m_atomic_species_MOD_rmaxsp
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Float64, (nspecies,))
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

function get_vcln()
    symbol = :__m_atomic_species_MOD_vcln
    nrspmax  = get_nrspmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (nrspmax,nspecies))
end

function get_vrsp()
    symbol = :__m_atomic_species_MOD_vrsp
    nrspmax  = get_nrspmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (nrspmax,nspecies))
end

# maximum allowed states for each species
function get_maxstsp()
    return 40 # hardcoded parameter
end

function get_evalsp()
    symbol = :__m_atomic_species_MOD_evalsp
    maxstsp = get_maxstsp()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (maxstsp,maxspecies))
end

function get_spzn()
    symbol = :__m_atomic_species_MOD_spzn
    nspecies = get_nspecies()
    maxspecies = get_maxspecies()
    # only return 1:nspecies
    return _load_automatic_array(symbol, Float64, (maxspecies,))[1:nspecies]
end
