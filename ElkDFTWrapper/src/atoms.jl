function get_natmtot()
    natmtot = unsafe_load(cglobal( (:__m_atoms_MOD_natmtot, LIBLAPW), Int32 )) |> Int64
    return natmtot
end

