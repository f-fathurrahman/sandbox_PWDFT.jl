#
# Muffin tins radial grid and angular momentum variables
#
function get_lmaxapw()
    lmaxapw = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxapw, LIBLAPW), Int32 )) |> Int64
    return lmaxapw
end

function get_lmmaxapw()
    lmmaxapw = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmmaxapw, LIBLAPW), Int32 )) |> Int64
    return lmmaxapw
end

function get_lmaxo()
    lmaxo = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxo, LIBLAPW), Int32 )) |> Int64
    return lmaxo
end

function get_lmmaxo()
    lmmaxo = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmmaxo, LIBLAPW), Int32 )) |> Int64
    return lmmaxo
end


function get_lmaxi()
    lmaxi = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxi, LIBLAPW), Int32 )) |> Int64
    return lmaxi
end

function get_lmmaxi()
    lmmaxi = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmmaxi, LIBLAPW), Int32 )) |> Int64
    return lmmaxi
end


function get_nrmtmax()
    nrmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_nrmtmax, LIBLAPW), Int32 )) |> Int64
    return nrmtmax
end

function get_npmtmax()
    npmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_npmtmax, LIBLAPW), Int32 )) |> Int64
    return npmtmax
end

function get_npcmtmax()
    npcmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_npcmtmax, LIBLAPW), Int32 )) |> Int64
    return npcmtmax
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