function get_efermi()
    return unsafe_load(cglobal( (:__m_states_MOD_efermi, LIBLAPW), Float64 ))
end

function get_swidth()
    return unsafe_load(cglobal( (:__m_states_MOD_swidth, LIBLAPW), Float64 ))
end
