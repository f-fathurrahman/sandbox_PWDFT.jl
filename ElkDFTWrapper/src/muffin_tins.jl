#
# Muffin tins radial grid and angular momentum variables
#
function get_lmaxo()
    lmaxo = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxo, LIBLAPW), Int32 )) |> Int64
    return lmaxo
end

function get_nrmtmax()
    nrmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_nrmtmax, LIBLAPW), Int32 )) |> Int64
    return nrmtmax
end

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


# FIXME: Test this !!!!
function get_nrmt()
    nspecies = get_nspecies()
    ptr = cglobal( (:__m_muffin_tins_MOD_nrmt, LIBLAPW), Ptr{Int32} )
    nrmt = zeros(Int64, nspecies)
    # XXX: Using two unsafe_load?
    for i in 1:nspecies
        nrmt[i] = Int64(unsafe_load(ptr,i))
    end
    return nrmt
end


function get_nrmti()
    nspecies = unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
    ptr = cglobal( (:__m_muffin_tins_MOD_nrmti, LIBLAPW), Ptr{Int32} )
    nrmti = zeros(Int64,nspecies)
    for i in 1:nspecies
        nrmti[i] = Int64(unsafe_load(ptr,i))
    end
    return nrmti
end
