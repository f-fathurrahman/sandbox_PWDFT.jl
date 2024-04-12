# We will use Int64 instead of Int32 here

function get_natmtot()
    return unsafe_load(cglobal((:__m_atoms_MOD_natmtot, LIBLAPW), Int32 )) |> Int64
end

function get_natmmax()
    return unsafe_load(cglobal((:__m_atoms_MOD_natmmax, LIBLAPW), Int32 )) |> Int64
end

function get_nspecies()
    return unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
end

# This is hardcoded
function get_maxatoms()
    return 200
end

# This is hardcoded
function get_maxspecies()
    return 8
end

function get_natoms()
    symbol = :__m_atoms_MOD_natoms
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Int64, (maxspecies,))
end

function get_atposl()
    symbol = :__m_atoms_MOD_atposl
    maxatoms = get_maxatoms()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (3,maxatoms,maxspecies))
end

function get_atposc()
    symbol = :__m_atoms_MOD_atposc
    maxatoms = get_maxatoms()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (3,maxatoms,maxspecies))
end
