# tolerance for error in total charge
function get_epschg()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_epschg, LIBLAPW), Float64))
end

# total nuclear charge
function get_chgzn()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgzn, LIBLAPW), Float64))
end

# total valence charge
function get_chgval()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgval, LIBLAPW), Float64))
end

# total charge
function get_chgtot()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgtot, LIBLAPW), Float64))
end

# total core charge
function get_chgcrtot()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgcrtot, LIBLAPW), Float64))
end

# excess charge
function get_chgexs()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgexs, LIBLAPW), Float64))
end

# calculated total charge
function get_chgcalc()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgcalc, LIBLAPW), Float64))
end

# interstitial region charge
function get_chgir()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgir, LIBLAPW), Float64))
end

# total muffin-tin charge
function get_chgmttot()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgmttot, LIBLAPW), Float64))
end

# effective Wigner radius
function get_rwigner()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_rwigner, LIBLAPW), Float64))
end

# total moment magnitude
function get_momtotm()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_momtotm, LIBLAPW), Float64))
end

# total current magnitude
function get_curtotm()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_curtotm, LIBLAPW), Float64))
end

# core charges
function get_chgcr()
    symbol = :__m_charge_moment_current_MOD_chgcr
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (maxspecies,))
end

# core leakage charge
function get_chgcrlk()
    symbol = :__m_charge_moment_current_MOD_chgcrlk
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (natmtot,))
end

# muffin-tin charges
function get_chgmt()
    symbol = :__m_charge_moment_current_MOD_chgmt
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (natmtot,))
end

# muffin-tin moments
function get_mommt()
    symbol = :__m_charge_moment_current_MOD_mommt
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (3,natmtot))
end

# interstitial region moment
function get_momir()
    symbol = :__m_charge_moment_current_MOD_momir
    return _load_automatic_array(symbol, Float64, (3,))
end

# total muffin-tin moment
function get_mommttot()
    symbol = :__m_charge_moment_current_MOD_mommttot
    return _load_automatic_array(symbol, Float64, (3,))
end

# total current
function get_curtot()
    symbol = :__m_charge_moment_current_MOD_curtot
    return _load_automatic_array(symbol, Float64, (3,))
end

# total moment
function  get_momtot()
    symbol = :__m_charge_moment_current_MOD_momtot
    return _load_automatic_array(symbol, Float64, (3,))
end


