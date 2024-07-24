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

# APW-APW Hamiltonian integrals
#real(8), allocatable :: haa(:,:,:,:,:,:)
function get_haa()
    #ALLOCATE( haa(lmmaxo,apwordmax,0:lmaxapw,apwordmax,0:lmaxapw,natmtot) )
    symbol = :__m_hamiltonian_MOD_haa
    lmmaxo = get_lmmaxo()
    apwordmax = get_apwordmax()
    lmaxapw = get_lmaxapw()
    natmtot = get_natmtot()    
    haa = _load_allocatable_array(symbol, Float64, (lmmaxo,apwordmax,lmaxapw+1,apwordmax,lmaxapw+1,natmtot))
    return OffsetArray(haa, 1:lmmaxo, 1:apwordmax, 0:lmaxapw, 1:apwordmax, 0:lmaxapw, 1:natmtot)
end


#=
# index to the position of the local-orbitals in the H and O matrices
integer, allocatable :: idxlo(:,:,:)

# APW-local-orbital overlap integrals
real(8), allocatable :: oalo(:,:,:)

# local-orbital-local-orbital overlap integrals
real(8), allocatable :: ololo(:,:,:)

# local-orbital-APW Hamiltonian integrals
real(8), allocatable :: hloa(:,:,:,:,:)

# local-orbital-local-orbital Hamiltonian integrals
real(8), allocatable :: hlolo(:,:,:,:)

# complex Gaunt coefficient array
complex(8), allocatable :: gntyry(:,:,:)

# tefvr is .true. if the first-variational eigenvalue equation is to be solved
# as a real symmetric problem
logical tefvr

# tefvit is .true. if the first-variational eigenvalue equation is to be solved
# iteratively
logical tefvit

# minimum and maximum allowed number of eigenvalue equation iterations
integer minitefv,maxitefv

# eigenvalue mixing parameter for iterative solver
real(8) befvit

# iterative solver convergence tolerance
real(8) epsefvit

# type of eigenvalue solver to be used
integer evtype
=#