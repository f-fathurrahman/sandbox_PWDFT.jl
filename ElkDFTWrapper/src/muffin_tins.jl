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
function get_rlmt()
    #
    lmaxo = get_lmaxo()
    nrmtmax = get_nrmtmax()
    nspecies = get_nspecies()
    #
    Ndim1 = nrmtmax
    Ndim2 = 2*(lmaxo+2)
    Ndim3 = nspecies
    rlmt = zeros(Float64, Ndim1*Ndim2*Ndim3)
    ptr_rlmt = cglobal( (:__m_muffin_tins_MOD_rlmt, LIBLAPW), Ptr{Float64} )
    ip = 1
    for k in 1:Ndim3, j in 1:Ndim2, i in 1:Ndim1
        rlmt[ip] = unsafe_load(unsafe_load(ptr_rlmt,1),ip)
        ip = ip + 1
    end
    rlmt = reshape(rlmt, (Ndim1, Ndim2, Ndim3))
    # Finally, convert to OffsetArray
    rlmt = OffsetArray(rlmt, 1:nrmtmax, -lmaxo-1:lmaxo+2, 1:nspecies)
    return rlmt
end


# nrmt is an automatic array of size maxspecies
# we don't load it as Ptr{Int32}, instead it is Int32
function get_nrmt()
    nspecies = get_nspecies()
    ptr = cglobal( (:__m_muffin_tins_MOD_nrmt, LIBLAPW), Int32 )
    nrmt = zeros(Int64, nspecies)
    for i in 1:nspecies
        nrmt[i] = unsafe_load(ptr,i) |> Int64
    end
    return nrmt
end

# nrmti is an automatic array of size maxspecies
function get_nrmti()
    nspecies = get_nspecies()
    ptr = cglobal( (:__m_muffin_tins_MOD_nrmti, LIBLAPW), Int32 )
    nrmti = zeros(Int64,nspecies)
    for i in 1:nspecies
        nrmti[i] = unsafe_load(ptr,i) |> Int64
    end
    return nrmti
end


# npmt is an automatic array of size maxspecies
function get_npmt()
    nspecies = get_nspecies()
    ptr = cglobal( (:__m_muffin_tins_MOD_npmt, LIBLAPW), Int32 )
    npmt = zeros(Int64,nspecies)
    for i in 1:nspecies
        npmt[i] = unsafe_load(ptr,i) |> Int64
    end
    return npmt
end


# npmt is an automatic array of size maxspecies
function get_npmti()
    nspecies = get_nspecies()
    ptr = cglobal( (:__m_muffin_tins_MOD_npmti, LIBLAPW), Int32 )
    npmti = zeros(Int64,nspecies)
    for i in 1:nspecies
        npmti[i] = unsafe_load(ptr,i) |> Int64
    end
    return npmti
end


# Coarse packed MT arrays
# npcmt is an automatic array of size maxspecies
function get_npcmt()
    nspecies = get_nspecies()
    ptr = cglobal( (:__m_muffin_tins_MOD_npcmt, LIBLAPW), Int32 )
    npcmt = zeros(Int64,nspecies)
    for i in 1:nspecies
        npcmt[i] = unsafe_load(ptr,i) |> Int64
    end
    return npcmt
end


# npcmti is an automatic array of size maxspecies
function get_npcmti()
    nspecies = get_nspecies()
    ptr = cglobal( (:__m_muffin_tins_MOD_npcmti, LIBLAPW), Int32 )
    npcmti = zeros(Int64,nspecies)
    for i in 1:nspecies
        npcmti[i] = unsafe_load(ptr,i) |> Int64
    end
    return npcmti
end