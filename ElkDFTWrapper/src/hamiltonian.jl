# overlap and Hamiltonian matrices sizes at each k-point
function get_nmat()
    symbol = :__m_hamiltonian_MOD_nmat
    nspnfv = get_nspnfv()
    nkpt = get_nkpt()
    return _load_allocatable_array(symbol, Int64, (nspnfv,nkpt))
end

# maximum nmat over all k-points
function get_nmatmax()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_nmatmax, LIBLAPW), Int32)) |> Int64
end

# tefvr is .true. if the first-variational eigenvalue equation is to be solved
# as a real symmetric problem
function get_tefvr()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_tefvr, LIBLAPW), Bool))
end

# tefvit is .true. if the first-variational eigenvalue equation is to be solved iteratively
function get_tefvit()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_tefvit, LIBLAPW), Bool))
end

# minimum and maximum allowed number of eigenvalue equation iterations
function get_minitefv()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_minitefv, LIBLAPW), Int32)) |> Int64
end

function get_maxitefv()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_maxitefv, LIBLAPW), Int32)) |> Int64
end

# eigenvalue mixing parameter for iterative solver
function get_befvit()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_befvit, LIBLAPW), Float64))
end


# iterative solver convergence tolerance
function get_epsefvit()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_epsefvit, LIBLAPW), Float64))
end

# type of eigenvalue solver to be used
function get_evtype()
    return unsafe_load(cglobal((:__m_hamiltonian_MOD_evtype, LIBLAPW), Int32)) |> Int64
end


# index to the position of the local-orbitals in the H and O matrices
function get_idxlo()
    # idxlo(lolmmax,nlomax,natmtot)
    symbol = :__m_hamiltonian_MOD_idxlo
    lolmmax = get_lolmmax()
    nlomax = get_nlomax()
    natmtot = get_natmtot()    
    return _load_allocatable_array(symbol, Int64, (lolmmax, nlomax, natmtot))
end


# APW-APW Hamiltonian integrals
function get_haa()
    # haa(lmmaxo,apwordmax,0:lmaxapw,apwordmax,0:lmaxapw,natmtot)
    symbol = :__m_hamiltonian_MOD_haa
    lmmaxo = get_lmmaxo()
    apwordmax = get_apwordmax()
    lmaxapw = get_lmaxapw()
    natmtot = get_natmtot()
    haa = _load_allocatable_array(symbol, Float64, (lmmaxo,apwordmax,lmaxapw+1,apwordmax,lmaxapw+1,natmtot))
    return OffsetArray(haa, 1:lmmaxo, 1:apwordmax, 0:lmaxapw, 1:apwordmax, 0:lmaxapw, 1:natmtot)
end

# local-orbital-local-orbital Hamiltonian integrals
function get_hlolo()
    # hlolo(lmmaxo,nlomax,nlomax,natmtot)
    symbol = :__m_hamiltonian_MOD_hlolo
    lmmaxo = get_lmmaxo()
    nlomax = get_nlomax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (lmmaxo, nlomax, nlomax, natmtot))
end

# local-orbital-local-orbital overlap integrals
function get_ololo()
    # ololo(nlomax,nlomax,natmtot)
    symbol = :__m_hamiltonian_MOD_ololo
    nlomax = get_nlomax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (nlomax, nlomax, natmtot))
end

# APW-local-orbital overlap integrals
function get_oalo()
    # oalo(apwordmax,nlomax,natmtot)
    symbol = :__m_hamiltonian_MOD_oalo
    apwordmax = get_apwordmax()
    nlomax = get_nlomax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (apwordmax, nlomax, natmtot))
end

# local-orbital-APW Hamiltonian integrals
function get_hloa()
    #hloa(lmmaxo,apwordmax,0:lmaxapw,nlomax,natmtot))
    symbol = :__m_hamiltonian_MOD_hloa
    lmmaxo = get_lmmaxo()
    apwordmax = get_apwordmax()
    lmaxapw = get_lmaxapw()
    nlomax = get_nlomax()
    natmtot = get_natmtot()
    hloa = _load_allocatable_array(symbol, Float64, (lmmaxo, apwordmax, lmaxapw+1, nlomax, natmtot))
    return OffsetArray(hloa, 1:lmmaxo, 1:apwordmax, 0:lmaxapw, 1:nlomax, 1:natmtot)
end

# complex Gaunt coefficient array
function get_gntyry()
    # gntyry(lmmaxo,lmmaxapw,lmmaxapw) 
    symbol = :__m_hamiltonian_MOD_gntyry
    lmmaxo = get_lmmaxo()
    lmmaxapw = get_lmmaxapw()
    return _load_allocatable_array(symbol, ComplexF64, (lmmaxo, lmmaxapw, lmmaxapw))
end
