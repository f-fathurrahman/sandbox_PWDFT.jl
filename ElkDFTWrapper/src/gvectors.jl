function get_ngtot()
    ngtot = unsafe_load(cglobal((:__m_gvectors_MOD_ngtot, LIBLAPW), Int32)) |> Int64
    return ngtot
end