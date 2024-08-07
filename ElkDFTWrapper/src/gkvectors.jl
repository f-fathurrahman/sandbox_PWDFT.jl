# species for which the muffin-tin radius will be used for calculating gkmax
function get_isgkmax()
    return unsafe_load(cglobal((:__m_gkvectors_MOD_isgkmax, LIBLAPW), Int32)) |> Int64
end

# maximum number of G+k-vectors over all k-points
function get_ngkmax()
    return unsafe_load(cglobal((:__m_gkvectors_MOD_ngkmax, LIBLAPW), Int32)) |> Int64
end

# smallest muffin-tin radius times gkmax
function get_rgkmax()
    return unsafe_load(cglobal((:__m_gkvectors_MOD_rgkmax, LIBLAPW), Float64))
end

# maximum |G+k| cut-off for APW functions
function get_gkmax()
    return unsafe_load(cglobal((:__m_gkvectors_MOD_gkmax, LIBLAPW), Float64))
end

# number of G+k-vectors for augmented plane waves
function get_ngk()
    symbol = :__m_gkvectors_MOD_ngk
    nkpt = get_nkpt() # nppt
    nspnfv = get_nspnfv()
    return _load_allocatable_array(symbol, Int64, (nspnfv,nkpt))
end

# index from G+k-vectors to G-vectors
# igkig(ngkmax,nspnfv,nppt)
function get_igkig()
    symbol = :__m_gkvectors_MOD_igkig
    ngkmax = get_ngkmax()
    nkpt = get_nkpt() # nppt
    nspnfv = get_nspnfv()
    return _load_allocatable_array(symbol, Int64, (ngkmax,nspnfv,nkpt))
end

# G+k-vectors in lattice coordinates
# vgkl(3,ngkmax,nspnfv,nppt)
function get_vgkl()
    symbol = :__m_gkvectors_MOD_vgkl
    ngkmax = get_ngkmax()
    nkpt = get_nkpt() # nppt
    nspnfv = get_nspnfv()
    return _load_allocatable_array(symbol, Float64, (3,ngkmax,nspnfv,nkpt))
end

# G+k-vectors in Cartesian coordinates
function get_vgkc()
    symbol = :__m_gkvectors_MOD_vgkc
    ngkmax = get_ngkmax()
    nkpt = get_nkpt() # nppt
    nspnfv = get_nspnfv()
    return _load_allocatable_array(symbol, Float64, (3,ngkmax,nspnfv,nkpt))
end

# length of G+k-vectors
function get_gkc()
    symbol = :__m_gkvectors_MOD_gkc
    ngkmax = get_ngkmax()
    nkpt = get_nkpt() # nppt
    nspnfv = get_nspnfv()
    return _load_allocatable_array(symbol, Float64, (ngkmax,nspnfv,nkpt))
end

# structure factors for the G+k-vectors
# sfacgk(ngkmax, natmtot, nspnfv,nppt)
function get_sfacgk()
    symbol = :__m_gkvectors_MOD_sfacgk
    ngkmax = get_ngkmax()
    natmtot = get_natmtot()
    nkpt = get_nkpt() # nppt
    nspnfv = get_nspnfv()
    return _load_allocatable_array(symbol, ComplexF64, (ngkmax,natmtot,nspnfv,nkpt))
end


#=
COMPLEX(8), ALLOCATABLE :: sfacgk(:,:,:,:)
=#