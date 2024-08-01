# Fermi energy for second-variational states (second variational ???)
function get_efermi()
    return unsafe_load(cglobal( (:__m_states_MOD_efermi, LIBLAPW), Float64 ))
end

# number of empty states per atom and spin
# XXX This is real(8) ?
function get_nempty0()
    return unsafe_load(cglobal( (:__m_states_MOD_nempty0, LIBLAPW), Float64 ))
end

# number of empty states
function get_nempty()
    return unsafe_load(cglobal( (:__m_states_MOD_nempty, LIBLAPW), Int32 )) |> Int64
end

# number of first-variational states
function get_nstfv()
    return unsafe_load(cglobal( (:__m_states_MOD_nstfv, LIBLAPW), Int32 )) |> Int64
end

# number of second-variational states
function get_nstsv()
    return unsafe_load(cglobal( (:__m_states_MOD_nstsv, LIBLAPW), Int32 )) |> Int64
end

# smearing type
function get_stype()
    return unsafe_load(cglobal( (:__m_states_MOD_stype, LIBLAPW), Int32 )) |> Int64
end

# smearing width
function get_swidth()
    return unsafe_load(cglobal( (:__m_states_MOD_swidth, LIBLAPW), Float64 ))
end

# autoswidth is .true. if the smearing width is to be determined automatically
function get_autoswidth()
    return unsafe_load(cglobal( (:__m_states_MOD_autoswidth, LIBLAPW), Bool ))
end

# effective mass used in smearing width formula
function get_mstar()
    return unsafe_load(cglobal( (:__m_states_MOD_mstar, LIBLAPW), Float64 ))
end

# maximum allowed occupancy (1 or 2)
function get_occmax()
    return unsafe_load(cglobal( (:__m_states_MOD_occmax, LIBLAPW), Float64 ))
end

# convergence tolerance for occupancies
function get_epsocc()
    return unsafe_load(cglobal( (:__m_states_MOD_epsocc, LIBLAPW), Float64 ))
end

# second-variational occupation numbers
function get_occsv()
    symbol = :__m_states_MOD_occsv
    nstsv = get_nstsv()
    nkpt = get_nkpt()
    return _load_allocatable_array(symbol, Float64, (nstsv,nkpt))
end

# scissor correction applied when computing response functions
function get_scissor()
    return unsafe_load(cglobal( (:__m_states_MOD_scissor, LIBLAPW), Float64 ))
end

# density of states at the Fermi energy
function get_fermidos()
    return unsafe_load(cglobal( (:__m_states_MOD_fermidos, LIBLAPW), Float64 ))
end

# estimated indirect and direct band gaps
function get_bandgap()
    symbol = :__m_states_MOD_bandgap
    return _load_automatic_array(symbol, Float64, (2,))
end

# k-points of indirect and direct gaps
function get_ikgap()
    symbol = :__m_states_MOD_ikgap
    return _load_automatic_array(symbol, Int64, (3,))
end

# error tolerance for the first-variational eigenvalues
function get_evaltol()
    return unsafe_load(cglobal( (:__m_states_MOD_evaltol, LIBLAPW), Float64 ))
end

# second-variational eigenvalues
function get_evalsv()
    symbol = :__m_states_MOD_evalsv
    nstsv = get_nstsv()
    nkpt = get_nkpt()
    return _load_allocatable_array(symbol, Float64, (nstsv,nkpt))
end

# tevecsv is .true. if second-variational eigenvectors are calculated
function get_tevecsv()
    return unsafe_load(cglobal( (:__m_states_MOD_tevecsv, LIBLAPW), Bool ))
end

# maximum number of k-point and states indices in user-defined list
function get_maxkst()
    return 20 # parameter
end

# number of k-point and states indices in user-defined list
function get_nkstlist()
    return unsafe_load(cglobal( (:__m_states_MOD_nkstlist, LIBLAPW), Int32 )) |> Int64
end

# user-defined list of k-point and state indices
function get_kstlist()
    symbol = :__m_states_MOD_nkstlist
    maxkst = get_maxkst()
    return _load_automatic_array(symbol, Int64, (2,maxkst))
end
