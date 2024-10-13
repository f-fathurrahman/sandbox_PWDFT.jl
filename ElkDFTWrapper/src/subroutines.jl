# Need to call this before doing anything else with Elk
function call_read_input()
    ccall( (:read_input_, LIBLAPW), Cvoid, () )
    return
end

# This is the original gndstate
function call_gndstate()
    ccall( (:gndstate_, LIBLAPW), Cvoid, () )
    return
end


function call_my_gndstate(maxscl)
    ccall( (:my_gndstate_, LIBLAPW), Cvoid, 
        (Ref{Int32}, ),
        Int32(maxscl)
    )
    return
end

function call_my_gndstate_setup_mixing()
    ccall( (:my_gndstate_setup_mixing_, LIBLAPW), Cvoid, () )
    return
end

function call_my_gndstate_do_mixing()
    ccall( (:my_gndstate_do_mixing_, LIBLAPW), Cvoid, () )
    return
end

function call_my_gndstate_increment_iscl()
    ccall( (:my_gndstate_increment_iscl_, LIBLAPW), Cvoid, () )
    return
end

function call_rotrfmt!(rot, nr, nri, rfmt1, rfmt2)
    #CALL my_rotrfmt(symlatc(:,:,lspl), nr(is), nri(is), rfmt1(:,ja), rfmt2)
    # SUBROUTINE my_rotrfmt(rot, nr, nri, rfmt1, rfmt2)
    ccall( (:my_rotrfmt_, LIBLAPW), Cvoid,
        (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        rot, Int32(nr), Int32(nri), rfmt1, rfmt2
    )
end

function call_symrf!(vxcmt, vxcir)
    nrmt = get_nrmt()
    nrmti = get_nrmti()
    npmt = get_npmt()
    npmtmax = get_npmtmax()
    # call the my_* version of symrf (debug)
    # CALL my_symrf(nrmt, nrmti, npmt, npmtmax, vxcmt_, vxcir_)
    ccall( (:my_symrf_, LIBLAPW), Cvoid,
        (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        Int32.(nrmt), Int32.(nrmti), Int32.(npmt), Int32(npmtmax), vxcmt, vxcir
    )
end

function call_genvsig()
    ccall( (:genvsig_, LIBLAPW), Cvoid, () )
    return
end

function call_linengy()
    ccall( (:linengy_, LIBLAPW), Cvoid, () )
    return
end

function call_genapwfr()
    ccall( (:genapwfr_, LIBLAPW), Cvoid, () )
    return
end

function call_genlofr()
    ccall( (:genlofr_, LIBLAPW), Cvoid, () )
    return
end

# Call the original genapwlofr
function call_genapwlofr()
    ccall( (:genapwlofr_, LIBLAPW), Cvoid, () )
    return
end

function call_genevfsv()
    ccall( (:genevfsv_, LIBLAPW), Cvoid, () )
    return
end

function call_occupy()
    ccall( (:occupy_, LIBLAPW), Cvoid, () )
    return
end

# Debug
function call_my_occupy()
    ccall( (:my_occupy_, LIBLAPW), Cvoid, () )
    return
end

function call_rhomag()
    ccall( (:rhomag_, LIBLAPW), Cvoid, () )
    return
end

function call_energy()
    ccall( (:energy_, LIBLAPW), Cvoid, () )
    return
end

function call_symmetry()
    ccall( (:symmetry_, LIBLAPW), Cvoid, () )
    return
end

function call_rhoinit()
    ccall( (:rhoinit_, LIBLAPW), Cvoid, () )
    return
end

function call_maginit()
    ccall( (:maginit_, LIBLAPW), Cvoid, () )
    return
end

# FIXME: pass txc?
function call_potks(; txc=true)
    ccall( (:potks_, LIBLAPW), Cvoid, (Ref{Bool},), txc )
    return
end

function call_potks_no_symm(; txc=true)
    # The version called here has my_ as prefix
    ccall( (:my_potks_no_symm_, LIBLAPW), Cvoid, (Ref{Bool},), true )
    return
end

function call_potcoul()
    ccall( (:potcoul_, LIBLAPW), Cvoid, () )
    return
end

# We call potxc_default here instead of potxc as it is simpler
function call_potxc()
    ccall( (:potxc_default_, LIBLAPW), Cvoid, () )
    return
end

function call_gencore()
    ccall( (:my_gencore_, LIBLAPW), Cvoid, () )
    return
end

# Compute hlolo, hloa, haa
function call_hmlrad()
    ccall( (:hmlrad_, LIBLAPW), Cvoid, () )
    return
end

# ololo, oalo
function call_olprad()
    ccall( (:olprad_, LIBLAPW), Cvoid, () )
    return
end

function call_info_apwlo()
    ccall( (:info_apwlo_, LIBLAPW), Cvoid, () )
    return
end

#=

call_genevfsv()

  ALLOCATE(evalfv(nstfv,nspnfv))
  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv), evecsv(nstsv,nstsv))
  # second variational states probably are not needed
  DO ik=1,nkpt
    ! solve the first- and second-variational eigenvalue equations
    CALL eveqn(ik, evalfv, evecfv, evecsv)
  ENDDO

call_eveqn  <-- for each k-points
SUBROUTINE eveqn(ik, evalfv, evecfv, evecsv)

evecsv <--- not needed ?

# matching coefs?
  ALLOCATE(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  
  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)

DO jspn=1,nspnfv

#  find the matching coefficients
    CALL match(ngk(jspn,ik),vgkc(:,:,jspn,ik),gkc(:,jspn,ik), &
     sfacgk(:,:,jspn,ik),apwalm(:,:,:,:,jspn))
    
# solve the first-variational eigenvalue equation

# not iterative

! directly
      CALL eveqnfv(
      nmat(jspn,ik), ngk(jspn,ik), igkig(:,jspn,ik), 
      vkc(:,ik), &
       vgkc(:,:,jspn,ik), apwalm(:,:,:,:,jspn), evalfv(:,jspn),
    evecfv(:,:,jspn))

  ENDDO 

=#


# TODO: make a better interface for this
function call_atom!(
    SOL, ptnucl,
    zn, nst_,
    n_, l_, k_, occ,
    xctype_, xcgrad,
    nr_, r,
    evals, rho, vr, rwf
)
# SOL: 

#  REAL(8), intent(in) :: sol speed of light (in atomic unit)
#  LOGICAL, intent(in) :: ptnucl
#  REAL(8), intent(in) :: zn
#  INTEGER, intent(in) :: nst
#  INTEGER, intent(in) :: n(nst),l(nst),k(nst)
#  REAL(8), intent(inout) :: occ(nst)
#  INTEGER, intent(in) :: xctype(3),xcgrad
#  INTEGER, intent(in) :: nr
#  REAL(8), intent(in) :: r(nr)
#  REAL(8), intent(out) :: eval(nst)
#  REAL(8), intent(out) :: rho(nr),vr(nr)
#  REAL(8), intent(out) :: rwf(nr,2,nst)

    # NOTE: Arguments with trailing underscore must be converted to Int32
    nst = convert(Int32, nst_)
    n = convert(Array{Int32}, n_)
    l = convert(Array{Int32}, l_)
    k = convert(Array{Int32}, k_)
    nr = convert(Int32, nr_)
    xctype = convert(Array{Int32}, xctype_)

    ccall(
        (:atom_, LIBLAPW),
        Cvoid, (
            Ref{Float64}, Ref{Bool},
            Ref{Float64}, Ref{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
            Ptr{Int32}, Ref{Bool},
            Ref{Int32}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        SOL, ptnucl,
        zn, nst,
        n, l, k, occ,
        xctype, xcgrad,
        nr, r,
        evals, rho, vr, rwf
    )

    return
end


function get_apwalm(ik)
    # COMPLEX(8) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
    ngkmax = get_ngkmax()
    apwordmax = get_apwordmax()
    lmmaxapw = get_lmmaxapw()
    natmtot = get_natmtot()
    nspnfv = get_nspnfv()
    apwalm = zeros(ComplexF64, ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv)
    # Loop over jspn is done in driver_match (Fortran)
    ccall( (:driver_match_, LIBLAPW), Cvoid,
        (Ref{Int32}, Ptr{ComplexF64}),
        Int32(ik), apwalm
    )
    return apwalm
end

function test_sbesseldm()

    # SUBROUTINE sbesseldm(m,lmax,x,djl)
    m = 2
    lmax = 4
    x = 0.1
    @info "m = $m lmax = $lmax x = $x"
    djl = zeros(Float64, lmax+1)
    ccall( (:sbesseldm_, LIBLAPW), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Float64}),
        Int32(m), Int32(lmax), x, djl
    )
    return djl
end


# This is utilized in hmlaa and olpaa
function call_zmctmu(tcr, a, b, c)
    l = size(a, 1)
    n = size(a, 2)
    @assert size(b, 1) == l
    @assert size(b, 2) == n
    ld = size(c, 1)
    ccall( (:zmctmu_, LIBLAPW), Cvoid,
        (Ref{Bool}, Ref{Int32}, Ref{Int32},
         Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Int32}, Ptr{ComplexF64}),
        tcr, Int32(l), Int32(n), a, b, Int32(ld), c 
    )
    return

#=
Typical values (for testing zmctmu) 
- tcr is tefvr - (e)igen (f)irst (v)ariational (r)eal
- l is lmoapw
- n is Ngk
- ld is nmatp, second dim of c is also nmatp

ik = 1
isp = 1
tefvr = elk.get_tefvr()
nmat = elk.get_nmat()
ngk = elk.get_ngk()
nmatk = nmat[ik]
lmoapw = elk.get_lmoapw()

n = ngk[ik]
nmatk = nmat[ik]
l = lmoapw[isp]

a = rand(ComplexF64, l, n)
b = rand(ComplexF64, l, n)
c = zeros(ComplexF64, nmatk, nmatk)
elk.call_zmctmu(tefvr, a, b, c)
=#

end


function get_hamiltonian(; ispin=1, ik=1)

    ngkmax = get_ngkmax()
    apwordmax = get_apwordmax()
    lmmaxapw = get_lmmaxapw()
    natmtot = get_natmtot()
    nspnfv = get_nspnfv()
    apwalm = zeros(ComplexF64, ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv)

    nmat = get_nmat()
    N = nmat[ispin,ik]
    Hmat = zeros(ComplexF64, N, N)
    Omat = zeros(ComplexF64, N, N)
    ccall( (:debug_hamiltonian_, LIBLAPW), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}),
        Int32(ispin), Int32(ik), apwalm, Hmat, Omat
    )
    serialize("apwalm_ispin_$(ispin)_ik_$(ik).dat", apwalm)
    serialize("Hmat_ispin_$(ispin)_ik_$(ik).dat", Hmat)
    serialize("Omat_ispin_$(ispin)_ik_$(ik).dat", Omat)
    return apwalm, Hmat, Omat
end
