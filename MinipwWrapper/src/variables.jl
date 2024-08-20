function get_lmaxkb()
    lmaxkb = unsafe_load(cglobal((:__uspp_param_MOD_lmaxkb, LIBMINIPW), Int32)) |> Int64
    return lmaxkb
end

