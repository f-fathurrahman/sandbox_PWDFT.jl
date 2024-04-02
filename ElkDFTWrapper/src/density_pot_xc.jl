# this assumes that rhomt is already allocated or initialized, by calling rhoinit
function get_rhomt()
    # Dimensions
    nrmtmax = unsafe_load(cglobal( (:__m_muffin_tins_MOD_npmtmax, LIBLAPW), Int32 )) |> Int64
    natmtot = unsafe_load(cglobal( (:__m_atoms_MOD_natmtot, LIBLAPW), Int32 )) |> Int64
    # Read the array
    # rhomt is an allocatable array in Fortran
    ptr = cglobal( (:__m_density_pot_xc_MOD_rhomt, LIBLAPW), Ptr{Float64} )
    Ndim1 = nrmtmax
    Ndim2 = natmtot
    tmp = zeros(Float64,Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    rhomt = reshape(tmp, Ndim1, Ndim2)
    return rhomt
end

# get rhoir