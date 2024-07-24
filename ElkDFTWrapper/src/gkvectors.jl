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


#=
INTEGER, ALLOCATABLE :: ngk(:,:)

# index from G+k-vectors to G-vectors
INTEGER, ALLOCATABLE :: igkig(:,:,:)

# G+k-vectors in lattice coordinates
REAL(8), ALLOCATABLE :: vgkl(:,:,:,:)

# G+k-vectors in Cartesian coordinates
REAL(8), ALLOCATABLE :: vgkc(:,:,:,:)

# length of G+k-vectors
REAL(8), ALLOCATABLE :: gkc(:,:,:)

# structure factors for the G+k-vectors
COMPLEX(8), ALLOCATABLE :: sfacgk(:,:,:,:)
=#