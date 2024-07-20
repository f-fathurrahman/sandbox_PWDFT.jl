# G-vector grid sizes for coarse grid (G < 2*gkmax)
function get_ngdc()
    symbol = :__m_gvectors_MOD_ngdc
    return _load_automatic_array(symbol, Int64, (3,))
end

# G-vector grid sizes
function get_ngridg()
    symbol = :__m_gvectors_MOD_ngridg
    return _load_automatic_array(symbol, Int64, (3,))
end

# G-vector cut-off for interstitial potential and density
function get_gmaxvr()
    return unsafe_load(cglobal((:__m_gvectors_MOD_gmaxvr, LIBLAPW), Float64))
end

# total number of G-vectors
function get_ngtot()
    return unsafe_load(cglobal((:__m_gvectors_MOD_ngtot, LIBLAPW), Int32)) |> Int64
end

# total number of G-vectors for coarse grid (G < 2*gkmax)
function get_ngtc()
    return unsafe_load(cglobal((:__m_gvectors_MOD_ngtc, LIBLAPW), Int32)) |> Int64
end


# integer grid intervals for each direction
function get_intgv()
    symbol = :__m_gvectors_MOD_intgv
    return _load_automatic_array(symbol, Int64, (2,3))
end

# number of G-vectors with G < gmaxvr
function get_ngvec()
    return unsafe_load(cglobal((:__m_gvectors_MOD_ngvec, LIBLAPW), Int32)) |> Int64
end

# number of G-vectors for coarse grid (G < 2*gkmax)
function get_ngvc()
    return unsafe_load(cglobal((:__m_gvectors_MOD_ngvc, LIBLAPW), Int32)) |> Int64
end

# G-vector integer coordinates (i1,i2,i3)
function get_ivg()
    symbol = :__m_gvectors_MOD_ivg
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Int64, (3,ngtot))
end

# map from (i1,i2,i3) to G-vector index
function get_ivgig()
    symbol = :__m_gvectors_MOD_ivgig
    intgv = get_intgv()
    idx1 = intgv[1,1]:intgv[2,1]
    idx2 = intgv[1,2]:intgv[2,2]
    idx3 = intgv[1,3]:intgv[2,3]
    ivgig = _load_allocatable_array(symbol, Int64, (length(idx1), length(idx2), length(idx3)))
    return OffsetArray(ivgig, (idx1,idx2,idx3))
end

# map from G-vector index to FFT array
function get_igfft()
    symbol = :__m_gvectors_MOD_igfft
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Int64, (ngtot,))
end

# map from G-vector index to FFT array for coarse grid (G < 2*gkmax)
#integer, allocatable :: igfc(:)
# ffr: XXX igfc is not used? Not yet used in Elk-6.3.2

# G-vectors in Cartesian coordinates
function get_vgc()
    symbol = :__m_gvectors_MOD_vgc
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (3,ngtot))
end

# length of G-vectors
function get_gc()
    symbol = :__m_gvectors_MOD_gc
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

# Coulomb Green's function in G-space = 4 pi / G^2
function get_gclg()
    symbol = :__m_gvectors_MOD_gclg
    ngvec = get_ngvec()
    return _load_allocatable_array(symbol, Float64, (ngvec,))
end

# spherical Bessel functions j_l(|G|R_mt)
function get_jlgrmt()
    symbol = :__m_gvectors_MOD_jlgrmt
    lnpsd = get_lnpsd()
    ngvec = get_ngvec()
    nspecies = get_nspecies()
    jlgrmt = _load_allocatable_array(symbol, Float64, (lnpsd+1, ngvec, nspecies))
    return OffsetArray(jlgrmt, (0:lnpsd, 1:ngvec, 1:nspecies))
end


# spherical harmonics of the G-vectors
function get_ylmg()
    symbol = :__m_gvectors_MOD_ylmg
    lmmaxo = get_lmmaxo()
    ngvec = get_ngvec()
    return _load_allocatable_array(symbol, ComplexF64, (lmmaxo, ngvec))
end

# structure factors for the G-vectors
function get_sfacg()
    symbol = :__m_gvectors_MOD_sfacg
    ngvec = get_ngvec()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, ComplexF64, (ngvec, natmtot))
end

# smooth step function form factors for all species and G-vectors
function get_ffacg()
    symbol = :__m_gvectors_MOD_ffacg
    ngtot = get_ngtot()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (ngtot, nspecies))
end

# characteristic function in G-space: 0 inside the muffin-tins and 1 outside
function get_cfunig()
    symbol = :__m_gvectors_MOD_cfunig
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, ComplexF64, (ngtot,))
end

# characteristic function in real-space: 0 inside the muffin-tins and 1 outside
function get_cfunir()
    symbol = :__m_gvectors_MOD_cfunir
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

#=
XXX Not yet used in Elk-6.3.2
# characteristic function in real-space for coarse grid (G < 2*gkmax)
real(8), allocatable :: cfrc(:)
=#