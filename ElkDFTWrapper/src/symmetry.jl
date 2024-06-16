# maximum of symmetries allowed
function get_maxsymcrys()
    # Simply return a constant (should be the same as in Elk)
    # This is not exported as symbol in so lib (?), so we use hardcoded value
    return 192
end

# type of symmetry allowed for the crystal
#  0 : only the identity element is used
#  1 : full symmetry group is used
#  2 : only symmorphic symmetries are allowed
function get_symtype()
    return unsafe_load(cglobal((:__m_symmetry_MOD_symtype, LIBLAPW), Int32 )) |> Int64
end

# number of Bravais lattice point group symmetries
function get_nsymlat()
    return unsafe_load(cglobal((:__m_symmetry_MOD_nsymlat, LIBLAPW), Int32 )) |> Int64
end

# number of crystal symmetries
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

# XXX: symlat matrices are transpose of PWDFT.SymmetryInfo (?)
function get_symlat()
    symbol = :__m_symmetry_MOD_symlat
    return _load_automatic_array(symbol, Int64, (3,3,48))
end

# lattice point group symmetries in Cartesian coordinates
function get_symlatc()
    symbol = :__m_symmetry_MOD_symlatc
    return _load_automatic_array(symbol, Float64, (3,3,48))
end

# determinants of lattice symmetry matrices (1 or -1)
function get_symlatd()
    symbol = :__m_symmetry_MOD_symlatd
    return _load_automatic_array(symbol, Int64, (48,))
end

# index to inverses of the lattice symmetries
function get_isymlat()
    symbol = :__m_symmetry_MOD_isymlat
    return _load_automatic_array(symbol, Int64, (48,))
end

# eqatoms(ia,ja,is) is true if atoms ia and ja are equivalent
function get_eqatoms()
    symbol = :__m_symmetry_MOD_eqatoms
    natmmax = get_natmmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Bool, (natmmax,natmmax,nspecies))
end

# equivalent atom index for each crystal symmetry
function get_ieqatom()
    symbol = :__m_symmetry_MOD_ieqatom
    natmmax = get_natmmax()
    nspecies = get_nspecies()
    maxsymcrys = get_maxsymcrys()
    return _load_allocatable_array(symbol, Int64, (natmmax,nspecies,maxsymcrys))
end

# number of site symmetries
function get_nsymsite()
    symbol = :__m_symmetry_MOD_nsymsite
    natmmax = get_natmtot()
    return _load_allocatable_array(symbol, Int64, (natmmax,))
end

# site symmetry spatial rotation element in lattice point group
function get_lsplsyms() 
    symbol = :__m_symmetry_MOD_lsplsyms
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Int64, (48,natmtot))
end

# site symmetry global spin rotation element in lattice point group
function get_lspnsyms() 
    symbol = :__m_symmetry_MOD_lspnsyms
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Int64, (48,natmtot))
end