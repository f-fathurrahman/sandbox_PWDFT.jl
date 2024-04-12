function get_avec()
    symbol = :__m_lattice_MOD_avec
    return _load_automatic_array(symbol, Float64, (3,3))
end

function get_ainv()
    symbol = :__m_lattice_MOD_ainv
    return _load_automatic_array(symbol, Float64, (3,3))
end

function get_bvec()
    symbol = :__m_lattice_MOD_bvec
    return _load_automatic_array(symbol, Float64, (3,3))
end

function get_binv()
    symbol = :__m_lattice_MOD_binv
    return _load_automatic_array(symbol, Float64, (3,3))
end
