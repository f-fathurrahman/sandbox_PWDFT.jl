function get_maxscl()
    return unsafe_load(cglobal((:__m_convergence_MOD_maxscl, LIBLAPW), Int32)) |> Int64
end

function get_iscl()
    return unsafe_load(cglobal((:__m_convergence_MOD_iscl, LIBLAPW), Int32)) |> Int64
end

function get_epspot()
    return unsafe_load(cglobal((:__m_convergence_MOD_epspot, LIBLAPW), Float64))
end

function get_epsengy()
    return unsafe_load(cglobal((:__m_convergence_MOD_epsengy, LIBLAPW), Float64))
end

function get_epsforce()
    return unsafe_load(cglobal((:__m_convergence_MOD_epsforce, LIBLAPW), Float64))
end

function get_epsstress()
    return unsafe_load(cglobal((:__m_convergence_MOD_epsstress, LIBLAPW), Float64))
end

