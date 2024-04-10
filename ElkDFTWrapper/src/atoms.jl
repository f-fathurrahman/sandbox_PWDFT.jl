function get_natmtot()
    natmtot = unsafe_load(cglobal( (:__m_atoms_MOD_natmtot, LIBLAPW), Int32 )) |> Int64
    return natmtot
end

# We will use Int64 instead of Int32 here
function get_nspecies()
    nspecies = unsafe_load(cglobal((:__m_atoms_MOD_nspecies, LIBLAPW), Int32)) |> Int64
    return nspecies
end

