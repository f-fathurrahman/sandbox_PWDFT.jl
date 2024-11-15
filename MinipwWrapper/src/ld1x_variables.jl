# hardcoded, from radial_grids module
function get_ld1x_ndmx()
    return 3500   # the maximum mesh size 
end

# hardcoded, from ld1_parameters
function get_ld1x_nwfx()
    return 38
end

# the number of wavefunctions
function get_ld1x_nwf()
    return unsafe_load(cglobal((:__ld1inc_MOD_nwf, LIBMINIPW), Int32)) |> Int64
end

# 1 (default) or 2 (if lsd=true)
function get_ld1x_nspin()
    return unsafe_load(cglobal((:__ld1inc_MOD_nspin, LIBMINIPW), Int32)) |> Int64
end

# the main quantum number
function get_ld1x_nn()
    symbol = :__ld1inc_MOD_nn
    nwfx = get_ld1x_nwfx()
    return _load_automatic_array(symbol, Int64, (nwfx,))
end

# the orbital angular momentum
function get_ld1x_ll()
    symbol = :__ld1inc_MOD_ll
    nwfx = get_ld1x_nwfx()
    return _load_automatic_array(symbol, Int64, (nwfx,))
end

# spin of the wfc. if(.not.lsd) all 1 (default)
function get_ld1x_isw()
    symbol = :__ld1inc_MOD_isw
    nwfx = get_ld1x_nwfx()
    return _load_automatic_array(symbol, Int64, (nwfx,))
end

# the total angular momentum
function get_ld1x_jj()
    symbol = :__ld1inc_MOD_jj
    nwfx = get_ld1x_nwfx()
    return _load_automatic_array(symbol, Float64, (nwfx,))
end

# the occupations of the all-electron atom
function get_ld1x_oc()
    symbol = :__ld1inc_MOD_oc
    nwfx = get_ld1x_nwfx()
    return _load_automatic_array(symbol, Float64, (nwfx,))
end

# the ionic charge (Float64)
function get_ld1x_zed()
    return unsafe_load(cglobal((:__ld1inc_MOD_zed, LIBMINIPW), Float64))
end

# the number of electrons
function get_ld1x_enne()
    return unsafe_load(cglobal((:__ld1inc_MOD_enne, LIBMINIPW), Float64))
end
