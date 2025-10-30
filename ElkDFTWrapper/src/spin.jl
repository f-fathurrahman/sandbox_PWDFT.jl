# number of spin-dependent first-variational functions per state
function get_nspnfv()
    return unsafe_load(cglobal((:__m_spin_MOD_nspnfv, LIBLAPW), Int32)) |> Int64
end

# spinpol is .true. for spin-polarised calculations
function get_spinpol()
    return unsafe_load(cglobal((:__m_spin_MOD_spinpol, LIBLAPW), Bool))
end

# spinorb is .true. for spin-orbit coupling
function get_spinorb()
    return unsafe_load(cglobal((:__m_spin_MOD_spinorb, LIBLAPW), Bool))
end

# dimension of magnetisation and magnetic vector fields (1 or 3)
function get_ndmag()
    return unsafe_load(cglobal((:__m_spin_MOD_ndmag, LIBLAPW), Int32)) |> Int64
end

# second-variational spinor dimension (1 or 2)
function get_nspinor()
    return unsafe_load(cglobal((:__m_spin_MOD_nspinor, LIBLAPW), Int32)) |> Int64
end

# scale factor of spin-orbit coupling term in Hamiltonian
function get_socscf()
    return unsafe_load(cglobal((:__m_spin_MOD_socscf, LIBLAPW), Float64))
end

# ncmag is .true. if the magnetisation is non-collinear, i.e. when ndmag = 3
function get_ncmag()
    return unsafe_load(cglobal((:__m_spin_MOD_ncmag, LIBLAPW), Bool))
end


# if cmagz is .true. then collinear magnetism along the z-axis is enforced
function get_cmagz()
    return unsafe_load(cglobal((:__m_spin_MOD_cmagz, LIBLAPW), Bool))
end

# spcpl is .true. if the up and down spins are coupled
function get_spcpl()
    return unsafe_load(cglobal((:__m_spin_MOD_spcpl, LIBLAPW), Bool))
end

# fixed spin moment type
#  0      : none
#  1 (-1) : total moment (direction)
#  2 (-2) : individual muffin-tin moments (direction)
#  3 (-3) : total and muffin-tin moments (direction)
function get_fsmtype()
    return unsafe_load(cglobal((:__m_spin_MOD_fsmtype, LIBLAPW), Int32)) |> Int64
end

# fixed total spin magnetic moment
function get_momfix()
    symbol = :__m_spin_MOD_momfix
    return _load_automatic_array(symbol, Float64, (3,))
end

# fixed spin moment global effective field in Cartesian coordinates
function get_bfsmc()
    symbol = :__m_spin_MOD_bfsmc
    return _load_automatic_array(symbol, Float64, (3,))
end

# muffin-tin fixed spin moments
function get_mommtfx()
    symbol = :__m_spin_MOD_mommtfix
    maxatoms = get_maxatoms()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (3,maxatoms,maxspecies))
end

# global external magnetic field in Cartesian coordinates
function get_bfieldc()
    symbol = :__m_spin_MOD_bfieldc
    return _load_automatic_array(symbol, Float64, (3,))
end

# initial field
function get_bfieldc0()
    symbol = :__m_spin_MOD_bfieldc0
    return _load_automatic_array(symbol, Float64, (3,))
end

# external magnetic field in each muffin-tin in Cartesian coordinates
function get_bfcmt()
    symbol = :__m_spin_MOD_bfcmt
    maxatoms = get_maxatoms()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (3,maxatoms,maxspecies))
end

# initial field as read in from input file
function get_bfcmt0()
    symbol = :__m_spin_MOD_bfcmt0
    maxatoms = get_maxatoms()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (3,maxatoms,maxspecies))
end

#=
# muffin-tin fixed spin moment effective fields in Cartesian coordinates
real(8), allocatable :: bfsmcmt(:,:)

# fixed spin moment field step size
real(8) taufsm

# magnitude of random vectors added to muffin-tin fields
real(8) rndbfcmt

# external magnetic fields are multiplied by reducebf after each s.c. loop
real(8) reducebf

# spinsprl is .true. if a spin-spiral is to be calculated
logical spinsprl

# ssdph is .true. if the muffin-tin spin-spiral magnetisation is de-phased
logical ssdph

# map from second- to first-variational spin index
integer jspnfv(2)

# spin-spiral q-vector in lattice coordinates
real(8) vqlss(3)

# spin-spiral q-vector in Cartesian coordinates
real(8) vqcss(3)

# current q-point in spin-spiral supercell calculation
integer iqss

# number of primitive unit cells in spin-spiral supercell
integer nscss

# number of fixed spin direction points on the sphere for finding the magnetic
# anisotropy energy (MAE)
integer npmae0,npmae

# (theta,phi) coordinates for each MAE direction
real(8), allocatable :: tpmae(:,:)
=#