function get_maxsymcrys()
    return 192
end

# type of symmetry allowed for the crystal
#  0 : only the identity element is used
#  1 : full symmetry group is used
#  2 : only symmorphic symmetries are allowed
function get_symtype()
    return unsafe_load(cglobal((:__m_symmetry_MOD_symtype, LIBLAPW), Int32 )) |> Int64
end

function get_nsymlat()
    return unsafe_load(cglobal((:__m_symmetry_MOD_nsymlat, LIBLAPW), Int32 )) |> Int64
end

function get_nsymcrys()
    return unsafe_load(cglobal((:__m_symmetry_MOD_nsymcrys, LIBLAPW), Int32 )) |> Int64
end

# Translation vector in lattice coordinate
function get_vtlsymc()
    symbol = :__m_symmetry_MOD_vtlsymc
    maxsymcrys = get_maxsymcrys()
    return _load_automatic_array(symbol, Float64, (3,maxsymcrys))
end

# Translation vector in Cartesian coordinate
function get_vtcsymc()
    symbol = :__m_symmetry_MOD_vtcsymc
    maxsymcrys = get_maxsymcrys()
    return _load_automatic_array(symbol, Float64, (3,maxsymcrys))
end

function get_symlat()
    symbol = :__m_symmetry_MOD_symlat
    return _load_automatic_array(symbol, Int64, (3,3,48))
end

function get_symlatc()
    symbol = :__m_symmetry_MOD_symlatc
    return _load_automatic_array(symbol, Float64, (3,3,48))
end

function get_symlatd()
    symbol = :__m_symmetry_MOD_symlatd
    return _load_automatic_array(symbol, Int64, (48,))
end

function get_isymlat()
    symbol = :__m_symmetry_MOD_isymlat
    return _load_automatic_array(symbol, Int64, (48,))
end

function get_eqatoms()
    symbol = :__m_symmetry_MOD_eqatoms
    natmmax = get_natmmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Bool, (natmmax,natmmax,nspecies))
end

function get_ieqatom()
    symbol = :__m_symmetry_MOD_ieqatom
    natmmax = get_natmmax()
    nspecies = get_nspecies()
    maxsymcrys = get_maxsymcrys()
    return _load_allocatable_array(symbol, Int64, (natmmax,nspecies,maxsymcrys))
end

function get_nsymsite()
    symbol = :__m_symmetry_MOD_nsymsite
    natmmax = get_natmtot()
    return _load_allocatable_array(symbol, Int64, (natmmax,))
end
