# this assumes that rhomt is already allocated or initialized, by calling rhoinit
function get_rhomt()
    # Dimensions
    nrmtmax = get_nrmtmax()
    natmtot = get_natmtot()
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