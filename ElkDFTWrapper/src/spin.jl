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

#=
# scale factor of spin-orbit coupling term in Hamiltonian
real(8) socscf

# ncmag is .true. if the magnetisation is non-collinear, i.e. when ndmag = 3
logical ncmag

# if cmagz is .true. then collinear magnetism along the z-axis is enforced
logical cmagz

# spcpl is .true. if the up and down spins are coupled
logical spcpl

# fixed spin moment type
#  0      : none
#  1 (-1) : total moment (direction)
#  2 (-2) : individual muffin-tin moments (direction)
#  3 (-3) : total and muffin-tin moments (direction)
integer fsmtype

# fixed total spin magnetic moment
real(8) momfix(3)

# fixed spin moment global effective field in Cartesian coordinates
real(8) bfsmc(3)

# muffin-tin fixed spin moments
real(8) mommtfix(3,maxatoms,maxspecies)

# muffin-tin fixed spin moment effective fields in Cartesian coordinates
real(8), allocatable :: bfsmcmt(:,:)

# fixed spin moment field step size
real(8) taufsm

# global external magnetic field in Cartesian coordinates
real(8) bfieldc(3)

# initial field
real(8) bfieldc0(3)

# external magnetic field in each muffin-tin in Cartesian coordinates
real(8) bfcmt(3,maxatoms,maxspecies)

# initial field as read in from input file
real(8) bfcmt0(3,maxatoms,maxspecies)

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