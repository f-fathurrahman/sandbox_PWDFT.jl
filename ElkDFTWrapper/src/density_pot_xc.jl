# this assumes that rhomt is already allocated or initialized, by calling rhoinit
function get_rhomt()
    # Dimensions
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    # Read the array
    # rhomt is an allocatable array in Fortran
    ptr = cglobal( (:__m_density_pot_xc_MOD_rhomt, LIBLAPW), Ptr{Float64} )
    Ndim1 = npmtmax
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

function get_rhoir()
    # Dimensions
    ngtot = get_ngtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_rhoir, LIBLAPW), Ptr{Float64} )
    Ndim1 = ngtot
    tmp = zeros(Float64,Ndim1)
    ip = 1
    for i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    rhoir = tmp
    return rhoir
end

#
# Potentials, muffin tins
#

function get_vclmt()
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_vclmt, LIBLAPW), Ptr{Float64} )
    Ndim1 = npmtmax
    Ndim2 = natmtot
    tmp = zeros(Float64,Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    vclmt = reshape(tmp, Ndim1, Ndim2)
    return vclmt
end

function get_exmt()
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_exmt, LIBLAPW), Ptr{Float64} )
    Ndim1 = npmtmax
    Ndim2 = natmtot
    tmp = zeros(Float64,Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    exmt = reshape(tmp, Ndim1, Ndim2)
    return exmt
end

function get_ecmt()
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_ecmt, LIBLAPW), Ptr{Float64} )
    Ndim1 = npmtmax
    Ndim2 = natmtot
    tmp = zeros(Float64,Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    ecmt = reshape(tmp, Ndim1, Ndim2)
    return ecmt
end


function get_vxcmt()
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_vxcmt, LIBLAPW), Ptr{Float64} )
    Ndim1 = npmtmax
    Ndim2 = natmtot
    tmp = zeros(Float64,Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    vxcmt = reshape(tmp, Ndim1, Ndim2)
    return vxcmt
end

function get_vsmt()
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_vsmt, LIBLAPW), Ptr{Float64} )
    Ndim1 = npmtmax
    Ndim2 = natmtot
    tmp = zeros(Float64,Ndim1*Ndim2)
    ip = 1
    for j in 1:Ndim2, i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    vsmt = reshape(tmp, Ndim1, Ndim2)
    return vsmt
end


#
# Potentials, interstitial
#

function get_vclir()
    # Dimensions
    ngtot = get_ngtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_vclir, LIBLAPW), Ptr{Float64} )
    Ndim1 = ngtot
    tmp = zeros(Float64,Ndim1)
    ip = 1
    for i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    vclir = tmp
    return vclir
end

function get_exir()
    # Dimensions
    ngtot = get_ngtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_exir, LIBLAPW), Ptr{Float64} )
    Ndim1 = ngtot
    tmp = zeros(Float64,Ndim1)
    ip = 1
    for i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    exir = tmp
    return exir
end

function get_ecir()
    # Dimensions
    ngtot = get_ngtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_ecir, LIBLAPW), Ptr{Float64} )
    Ndim1 = ngtot
    tmp = zeros(Float64,Ndim1)
    ip = 1
    for i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    ecir = tmp
    return ecir
end


function get_vxcir()
    # Dimensions
    ngtot = get_ngtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_vxcir, LIBLAPW), Ptr{Float64} )
    Ndim1 = ngtot
    tmp = zeros(Float64,Ndim1)
    ip = 1
    for i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    vxcir = tmp
    return vxcir
end

function get_vsir()
    # Dimensions
    ngtot = get_ngtot()
    ptr = cglobal( (:__m_density_pot_xc_MOD_vsir, LIBLAPW), Ptr{Float64} )
    Ndim1 = ngtot
    tmp = zeros(Float64,Ndim1)
    ip = 1
    for i in 1:Ndim1
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    vsir = tmp
    return vsir
end
