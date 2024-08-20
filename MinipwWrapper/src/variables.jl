function get_lmaxkb()
    lmaxkb = unsafe_load(cglobal((:__uspp_param_MOD_lmaxkb, LIBMINIPW), Int32)) |> Int64
    return lmaxkb
end

# max non local angular momentum (l=0 to lmaxx)
function get_lmaxx()
    return 3 # HARCODED
end

# max number of angular momenta of Q
function get_lqmax()
    lmaxx = get_lmaxx()
    return 2*lmaxx + 1
end

function get_base_nsp()
    return unsafe_load(cglobal((:__ions_base_MOD_nsp, LIBMINIPW), Int32)) |> Int64
end
const get_Nspecies = get_base_nsp

function get_uspp_param_nh()
    symbol = :__uspp_param_MOD_nh
    nsp = get_base_nsp()
    return _load_automatic_array(symbol, Int64, (nsp,))
    # Int64 will be read as Int32, then converted to Int64 again
end

