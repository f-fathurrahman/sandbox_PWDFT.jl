function get_spincore()
    return unsafe_load(cglobal((:__m_core_states_MOD_spincore, LIBLAPW), Bool ))
    # Use Int32 ?
end

function get_nspncr()
    return unsafe_load(cglobal((:__m_core_states_MOD_nspncr, LIBLAPW), Int32 )) |> Int64
end

# occupancies for core states
function get_occcr()
    symbol = :__m_core_states_MOD_occcr
    nstspmax  = get_nstspmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (nstspmax,natmtot))
end

# eigenvalues for core states
function get_evalcr()
    symbol = :__m_core_states_MOD_evalcr
    nstspmax  = get_nstspmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (nstspmax,natmtot))
end


# radial wavefunctions for core states
function get_rwfcr()
    symbol = :__m_core_states_MOD_rwfcr
    nrspmax = get_nrspmax()
    nstspmax  = get_nstspmax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (nrspmax, 2, nstspmax, natmtot)) 
end


# radial charge density for core states
function get_rhocr()
    symbol = :__m_core_states_MOD_rhocr
    nrmtmax = get_nrmtmax()
    natmtot = get_natmtot()
    nspncr = get_nspncr()
    return _load_allocatable_array(symbol, Float64, (nrmtmax, natmtot, nspncr)) 
end

