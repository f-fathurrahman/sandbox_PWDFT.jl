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

function get_lsda_nspin()
    return unsafe_load(cglobal((:__lsda_mod_MOD_nspin, LIBMINIPW), Int32)) |> Int64
end

function get_dfftp_nnr()
    ccall((:__exposed_MOD_expose_dfftp_nnr, LIBMINIPW), Cvoid, ())
    return unsafe_load(cglobal((:__exposed_MOD_dfftp_nnr, LIBMINIPW), Int32)) |> Int64
end

function get_rho_of_r()
    symbol = :__exposed_MOD_rho_of_r
    nnr = get_dfftp_nnr()
    nspin = get_lsda_nspin()
    ccall((:__exposed_MOD_expose_rho_of_r, LIBMINIPW), Cvoid, ())
    return _load_allocatable_array(symbol, Float64, (nnr,nspin))
end
