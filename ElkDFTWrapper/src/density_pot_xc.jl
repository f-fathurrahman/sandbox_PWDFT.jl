# tssxc is .true. if scaled spin exchange-correlation (SSXC) is to be used
function get_tssxc()
    return unsafe_load(cglobal((:__m_density_pot_xc_MOD_tssxc, LIBLAPW), Bool))
end

# SSXC scaling factor
function get_ssxc()
    return unsafe_load(cglobal((:__m_density_pot_xc_MOD_ssxc, LIBLAPW), Float64))
end

# constant part of exchange-correlation potential
function get_vxc0()
    return unsafe_load(cglobal((:__m_density_pot_xc_MOD_vxc0, LIBLAPW), Float64))
end

# tbdip is .true. if the spin and current dipole fields are to be added to the
# Kohn-Sham magnetic field
function get_tbdip()
    return unsafe_load(cglobal((:__m_density_pot_xc_MOD_tbdip, LIBLAPW), Bool))
end

# trhonorm is .true. if the density is to be normalised after every iteration
# LOGICAL trhonorm
function get_trhonorm()
    return unsafe_load(cglobal((:__m_density_pot_xc_MOD_trhonorm, LIBLAPW), Bool))
end


function get_msmooth()
    return unsafe_load(cglobal((:__m_density_pot_xc_MOD_msmooth, LIBLAPW), Int32)) |> Int64
end

function get_lnpsd()
    return unsafe_load(cglobal((:__m_density_pot_xc_MOD_lnpsd, LIBLAPW), Int32)) |> Int64
end

function get_xctype()
    symbol = :__m_density_pot_xc_MOD_xctype
    return _load_automatic_array(symbol, Int64, (3,))    
end

# this assumes that rhomt is already allocated or initialized, by calling rhoinit
function get_rhomt()
    symbol = :__m_density_pot_xc_MOD_rhomt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot))
end

function get_rhoir()
    symbol = :__m_density_pot_xc_MOD_rhoir
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

#
# Potentials, muffin tins
#
function get_vclmt()
    symbol = :__m_density_pot_xc_MOD_vclmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot))
end

function get_exmt()
    symbol = :__m_density_pot_xc_MOD_exmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot))
end

function get_ecmt()
    symbol = :__m_density_pot_xc_MOD_ecmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot))
end


function get_vxcmt()
    symbol = :__m_density_pot_xc_MOD_vxcmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot))
end

function get_vsmt()
    symbol = :__m_density_pot_xc_MOD_vsmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot))
end


#
# Potentials, interstitial
#

function get_vclir()
    symbol = :__m_density_pot_xc_MOD_vclir
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

function get_exir()
    symbol = :__m_density_pot_xc_MOD_exir
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

function get_ecir()
    symbol = :__m_density_pot_xc_MOD_ecir
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

function get_vxcir()
    symbol = :__m_density_pot_xc_MOD_vxcir
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

function get_vsir()
    symbol = :__m_density_pot_xc_MOD_vsir
    ngtot = get_ngtot()
    return _load_allocatable_array(symbol, Float64, (ngtot,))
end

function get_vsig()
    symbol = :__m_density_pot_xc_MOD_vsig
    ngvec = get_ngvec()
    return _load_allocatable_array(symbol, ComplexF64, (ngvec,))
end



#
# muffin-tin and interstitial magnetisation vector field
function get_magmt()
    symbol = :__m_density_pot_xc_MOD_magmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ndmag = get_ndmag()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot,ndmag))
end

function get_magir()
    symbol = :__m_density_pot_xc_MOD_magir
    ngtot = get_ngtot()
    ndmag = get_ndmag()
    return _load_allocatable_array(symbol, Float64, (ngtot,ndmag))
end

# muffin-tin and interstitial exchange-correlation magnetic field
function get_bxcmt()
    symbol = :__m_density_pot_xc_MOD_bxcmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ndmag = get_ndmag()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot,ndmag))
end

function get_bxcir()
    symbol = :__m_density_pot_xc_MOD_bxcir
    ngtot = get_ngtot()
    ndmag = get_ndmag()
    return _load_allocatable_array(symbol, Float64, (ngtot,ndmag))
end

# muffin-tin Kohn-Sham effective magnetic field in spherical coordinates and on
# a coarse radial mesh
function get_bsmt()
    symbol = :__m_density_pot_xc_MOD_bsmt
    npmtmax = get_npmtmax()
    natmtot = get_natmtot()
    ndmag = get_ndmag()
    return _load_allocatable_array(symbol, Float64, (npmtmax,natmtot,ndmag))
end

function get_bsir()
    symbol = :__m_density_pot_xc_MOD_bsir
    ngtot = get_ngtot()
    ndmag = get_ndmag()
    return _load_allocatable_array(symbol, Float64, (ngtot,ndmag))
end