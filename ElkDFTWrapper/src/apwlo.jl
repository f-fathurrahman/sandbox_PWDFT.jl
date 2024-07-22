# energy step used for numerical calculation of energy derivatives
function get_deapwlo()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_deapwlo, LIBLAPW), Float64 ))
end

# maximum allowable APW order
function get_maxapword()
    return 4 # parameter, hardcoded
end

# maximum of apword over all angular momenta and species
function get_apwordmax()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_apwordmax, LIBLAPW), Int32 )) |> Int64
end

# polynomial order used for APW radial derivatives
function get_npapw()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_npapw, LIBLAPW), Int32 )) |> Int64
end

# maximum nlorb over all species
function get_nlomax()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_nlomax, LIBLAPW), Int32 )) |> Int64
end

# total number of local-orbitals
function get_nlotot()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_nlotot, LIBLAPW), Int32 )) |> Int64
end

# maximum number of local-orbitals
function get_maxlorb()
    return 200 # parameter, hardcoded
end

# maximum allowable local-orbital order
function get_maxlorbord()
    return 5 # parameter, hardcoded
end

#=

# band energy search tolerance
REAL(8) epsband

# maximum allowed change in energy during band energy search; enforced only if
# default energy is less than zero
REAL(8) demaxbnd

# minimum default linearisation energy over all APWs and local-orbitals
REAL(8) e0min

# maximum lorbord over all species
INTEGER :: lorbordmax

# polynomial order used for local-orbital radial derivatives
INTEGER :: nplorb

# maximum lorbl over all species
INTEGER lolmax

# (lolmax+1)^2
INTEGER lolmmax

# if autolinengy is .true. THEN  the fixed linearisation energies are set to the
# Fermi energy minus dlefe
LOGICAL autolinengy

# difference between linearisation and Fermi energies when autolinengy is .true.
REAL(8) dlefe

# lorbcnd is .true. if conduction state local-orbitals should be added
LOGICAL lorbcnd

# conduction state local-orbital order
INTEGER lorbordc

# excess order of the APW and local-orbital functions
INTEGER nxoapwlo

# excess local orbitals
INTEGER nxlo



# total number of APW coefficients (l, m and order) for each species
INTEGER :: lmoapw(maxspecies)

# APW order
INTEGER :: apword(0:maxlapw,maxspecies)

# APW initial linearisation energies
REAL(8) :: apwe0(maxapword,0:maxlapw,maxspecies)

# APW linearisation energies
REAL(8), ALLOCATABLE :: apwe(:,:,:)

# APW derivative order
INTEGER :: apwdm(maxapword,0:maxlapw,maxspecies)

# apwve is .true. if the linearisation energies are allowed to vary
LOGICAL :: apwve(maxapword,0:maxlapw,maxspecies)

# APW radial functions
REAL(8), ALLOCATABLE :: apwfr(:,:,:,:,:)

# derivate of radial functions at the muffin-tin surface
REAL(8), ALLOCATABLE :: apwdfr(:,:,:)

# number of local-orbitals
INTEGER :: nlorb(maxspecies)

# local-orbital order
INTEGER :: lorbord(maxlorb,maxspecies)

# local-orbital angular momentum
INTEGER :: lorbl(maxlorb,maxspecies)

# local-orbital initial energies
REAL(8) lorbe0(maxlorbord,maxlorb,maxspecies)

# local-orbital energies
REAL(8), ALLOCATABLE :: lorbe(:,:,:)

# local-orbital derivative order
INTEGER lorbdm(maxlorbord,maxlorb,maxspecies)

# lorbve is .true. if the linearisation energies are allowed to vary
LOGICAL lorbve(maxlorbord,maxlorb,maxspecies)

# local-orbital radial functions
REAL(8), ALLOCATABLE :: lofr(:,:,:,:)
=#

