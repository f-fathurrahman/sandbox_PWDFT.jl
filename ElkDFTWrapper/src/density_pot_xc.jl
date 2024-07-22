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
