#
# Muffin tins radial grid and angular momentum variables
#

function elk_write_rlmt()
    #
    lmaxo = unsafe_load(cglobal( (:__m_muffin_tins_MOD_lmaxo, LIBLAPW), Int32 )) |> Int64
    nrmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_nrmtmax, LIBLAPW), Int32 )) |> Int64
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
    rlmt = OffsetArray(rlmt, 1:nrmtmax, -lmaxo-1:lmaxo+2, 1:nspecies)
    
    serialize(joinpath(ELK_DATA_DIR, "rlmt.dat"), rlmt)
    return
end

