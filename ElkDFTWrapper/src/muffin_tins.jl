#
# Muffin tins radial grid and angular momentum variables
#

# maximum allowable angular momentum for augmented plane waves
function get_maxlapw()
    return 50 # hardcoded
end

# fraction of muffin-tin radius which constitutes the inner part
function get_fracinr()
    return unsafe_load(cglobal( (:__m_muffin_tins_MOD_fracinr, LIBLAPW), Float64 ))
end

# maximum angular momentum for augmented plane waves
# XXX: difference with maxlapw?
# This is the actual lmax used in APW basis expansion
function get_lmaxapw()
    lmaxapw = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxapw, LIBLAPW), Int32 )) |> Int64
    return lmaxapw
end

# (lmaxapw+1)^2
function get_lmmaxapw()
    lmmaxapw = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmmaxapw, LIBLAPW), Int32 )) |> Int64
    return lmmaxapw
end

# maximum angular momentum on the outer part of the muffin-tin
function get_lmaxo()
    lmaxo = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxo, LIBLAPW), Int32 )) |> Int64
    return lmaxo
end

# (lmaxo+1)^2
function get_lmmaxo()
    lmmaxo = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmmaxo, LIBLAPW), Int32 )) |> Int64
    return lmmaxo
end

# maximum angular momentum on the inner part of the muffin-tin
 function get_lmaxi()
    lmaxi = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxi, LIBLAPW), Int32 )) |> Int64
    return lmaxi
end

# (lmaxi+1)^2
function get_lmmaxi()
    lmmaxi = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmmaxi, LIBLAPW), Int32 )) |> Int64
    return lmmaxi
end

# maximum nrmt over all the species
function get_nrmtmax()
    nrmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_nrmtmax, LIBLAPW), Int32 )) |> Int64
    return nrmtmax
end

# maximum number of points over all packed muffin-tins (fine)
function get_npmtmax()
    npmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_npmtmax, LIBLAPW), Int32 )) |> Int64
    return npmtmax
end

# maximum number of points over all packed muffin-tins (coarse)
function get_npcmtmax()
    npcmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_npcmtmax, LIBLAPW), Int32 )) |> Int64
    return npcmtmax
end

function get_idxlm()
    symbol = :__m_muffin_tins_MOD_idxlm
    lmaxapw = get_lmaxapw()
    idxlm = _load_allocatable_array(symbol, Int64, (lmaxapw+1, 2*lmaxapw+1))
    return OffsetArray(idxlm, 0:lmaxapw, -lmaxapw:lmaxapw)
end

# lm -> l indices
function get_idxil()
    symbol = :__m_muffin_tins_MOD_idxil
    lmmaxapw = get_lmmaxapw()
    return _load_allocatable_array(symbol, Int64, (lmmaxapw,))
end

# lm -> m indices
function get_idxim()
    symbol = :__m_muffin_tins_MOD_idxim
    lmmaxapw = get_lmmaxapw()
    return _load_allocatable_array(symbol, Int64, (lmmaxapw,))
end


# XXX: Probably use a macro to generate these functions


# r^l on fine radial mesh (allocatable)
# This will be an OffsetArray
function get_rlmt()
    symbol = :__m_muffin_tins_MOD_rlmt
    #
    lmaxo = get_lmaxo()
    nrmtmax = get_nrmtmax()
    nspecies = get_nspecies()
    #
    Ndim1 = nrmtmax
    Ndim2 = 2*(lmaxo+2)
    Ndim3 = nspecies
    rlmt = _load_allocatable_array(symbol, Float64, (Ndim1,Ndim2,Ndim3))
    # Finally, convert to OffsetArray
    rlmt = OffsetArray(rlmt, 1:nrmtmax, -lmaxo-1:lmaxo+2, 1:nspecies)
    return rlmt
end

function get_rmt()
    symbol = :__m_muffin_tins_MOD_rmt
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Float64, (nspecies,))
end


# nrmt is an automatic array of size maxspecies
# we don't load it as Ptr{Int32}, instead it is Int32
function get_nrmt()
    symbol = :__m_muffin_tins_MOD_nrmt
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end

# nrmti is an automatic array of size maxspecies
function get_nrmti()
    symbol = :__m_muffin_tins_MOD_nrmti
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end

function get_nrcmti()
    symbol = :__m_muffin_tins_MOD_nrcmti
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end


# npmt is an automatic array of size maxspecies
function get_npmt()
    symbol = :__m_muffin_tins_MOD_npmt
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end


# npmt is an automatic array of size maxspecies
function get_npmti()
    symbol = :__m_muffin_tins_MOD_npmti
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end


# Coarse packed MT arrays
# npcmt is an automatic array of size maxspecies
function get_npcmt()
    symbol = :__m_muffin_tins_MOD_npcmt
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end


# npcmti is an automatic array of size maxspecies
function get_npcmti()
    symbol = :__m_muffin_tins_MOD_npcmti
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
end

# total muffin-tin volume
function get_omegamt()
    return unsafe_load(cglobal( (:__m_muffin_tins_MOD_omegamt, LIBLAPW), Float64 ))
end

# radial step length for coarse mesh
function get_lradstp()
    return unsafe_load(cglobal( (:__m_muffin_tins_MOD_lradstp, LIBLAPW), Int32 )) |> Int64
end

function get_nrcmtmax()
    return unsafe_load(cglobal( (:__m_muffin_tins_MOD_nrcmtmax, LIBLAPW), Int32 )) |> Int64
end

# number of coarse radial mesh points
function get_nrcmt()
    symbol = :__m_muffin_tins_MOD_nrcmt
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
    # XXX We only load nspecies data
end

# coarse muffin-tin radial mesh
function get_rcmt()
    symbol = :__m_muffin_tins_MOD_rcmt
    nspecies = get_nspecies()
    nrcmtmax = get_nrcmtmax()
    return _load_automatic_array(symbol, Float64, (nrcmtmax, nspecies))
end

# r^l on coarse radial mesh
function get_rlcmt()
    symbol = :__m_muffin_tins_MOD_rlcmt
    #
    lmaxo = get_lmaxo()
    nrcmtmax = get_nrcmtmax()
    nspecies = get_nspecies()
    #
    Ndim1 = nrcmtmax
    Ndim2 = 2*(lmaxo+2)
    Ndim3 = nspecies
    rlcmt = _load_allocatable_array(symbol, Float64, (Ndim1,Ndim2,Ndim3))
    return OffsetArray(rlcmt, 1:nrcmtmax, -lmaxo-1:lmaxo+2, 1:nspecies)
end


# weights for spline partial integration on fine radial mesh
function get_wprmt()
    symbol = :__m_muffin_tins_MOD_wprmt
    nrmtmax = get_nrmtmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (4,nrmtmax,nspecies))
end

# weights for spline partial integration on coarse radial mesh
function get_wprcmt()
    symbol = :__m_muffin_tins_MOD_wprcmt
    nrcmtmax = get_nrcmtmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (4,nrcmtmax,nspecies))
end

# weights for spline integration on coarse radial mesh
function get_wrcmt()
    symbol = :__m_muffin_tins_MOD_wrcmt
    nrcmtmax = get_nrcmtmax()
    nspecies = get_nspecies()
    return _load_allocatable_array(symbol, Float64, (nrcmtmax,nspecies))
end
