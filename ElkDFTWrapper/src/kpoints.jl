# total number of reduced k-points
function get_nkpt()
    return unsafe_load(cglobal((:__m_kpoints_MOD_nkpt, LIBLAPW), Int32)) |> Int64
end

# autokpt is .true. if the k-point set is determined automatically
function get_autokpt()
    return unsafe_load(cglobal((:__m_kpoints_MOD_autokpt, LIBLAPW), Bool))
end

# radius of sphere used to determine k-point density when autokpt is .true.
function get_radkpt()
    return unsafe_load(cglobal((:__m_kpoints_MOD_radkpt, LIBLAPW), Float64))
end

# type of reduction to perform on k-point set
#  0 : no reduction
#  1 : reduce with full crystal symmetry group
#  2 : reduce with symmorphic symmetries only
function get_reducek()
    return unsafe_load(cglobal((:__m_kpoints_MOD_reducek, LIBLAPW), Int32)) |> Int64
end

# total number of non-reduced k-points
function get_nkptnr()
    return unsafe_load(cglobal((:__m_kpoints_MOD_nkptnr, LIBLAPW), Int32)) |> Int64
end

# weight of each non-reduced k-point
function get_wkptnr()
    return unsafe_load(cglobal((:__m_kpoints_MOD_wkptnr, LIBLAPW), Float64))
end

# number of point group symmetries used for k-point reduction
function get_nsymkpt()
    return unsafe_load(cglobal((:__m_kpoints_MOD_nsymkpt, LIBLAPW), Int32)) |> Int64
end

# k-point grid sizes
#integer ngridk(3)
function get_ngridk()
    symbol = :__m_kpoints_MOD_ngridk
    return _load_automatic_array(symbol, Int64, (3,))
end

# k-point offset
#real(8) vkloff(3)
function get_vkloff()
    symbol = :__m_kpoints_MOD_vkloff
    return _load_automatic_array(symbol, Float64, (3,))
end

# corners of box in lattice coordinates containing the k-points, the zeroth
# vector is the origin
#real(8) kptboxl(3,0:3)
function get_kptboxl()
    symbol = :__m_kpoints_MOD_vkloff
    kptboxl = _load_automatic_array(symbol, Float64, (3,4))
    return OffsetArray(kptboxl, 1:3, 0:3)
end

# point group symmetry matrices used for k-point reduction
#integer symkpt(3,3,48)
function get_symkpt()
    symbol = :__m_kpoints_MOD_symkpt
    return _load_automatic_array(symbol, Int64, (3,3,48))
end

# locations of k-points on integer grid
#integer, allocatable :: ivk(:,:)
function get_ivk()
    symbol = :__m_kpoints_MOD_ivk
    nkptnr = get_nkptnr()
    return _load_allocatable_array(symbol, Int64, (3,nkptnr))
end

# map from (i1,i2,i3) to reduced k-point index
function get_ivkik()
    # ivkik(0:ngridk(1)-1, 0:ngridk(2)-1, 0:ngridk(3)-1)
    symbol = :__m_kpoints_MOD_ivkik
    ngridk = get_ngridk()
    ivkik = _load_allocatable_array(symbol, Int64, (ngridk[1], ngridk[2], ngridk[3]))
    return OffsetArray(ivkik, 0:ngridk[1]-1, 0:ngridk[2]-1, 0:ngridk[3]-1)
end

# map from (i1,i2,i3) to non-reduced k-point index
function get_ivkiknr()
    symbol = :__m_kpoints_MOD_ivkiknr
    ngridk = get_ngridk()
    ivkiknr = _load_allocatable_array(symbol, Int64, (ngridk[1], ngridk[2], ngridk[3]))
    return OffsetArray(ivkiknr, 0:ngridk[1]-1, 0:ngridk[2]-1, 0:ngridk[3]-1)
end

# k-points in lattice coordinates
function get_vkl()
    symbol = :__m_kpoints_MOD_vkl
    nkptnr = get_nkptnr()
    return _load_allocatable_array(symbol, Float64, (3,nkptnr))
end

# k-points in Cartesian coordinates
function get_vkc()
    symbol = :__m_kpoints_MOD_vkc
    nkptnr = get_nkptnr()
    return _load_allocatable_array(symbol, Float64, (3,nkptnr))
end

# reduced k-point weights
function get_wkpt()
    symbol = :__m_kpoints_MOD_wkpt
    nkptnr = get_nkptnr() # XXX the size is overallocated ?
    return _load_allocatable_array(symbol, Float64, (nkptnr,))
end


