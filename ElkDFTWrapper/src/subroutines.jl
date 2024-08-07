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
    ccall( (:my_linengy_, LIBLAPW), Cvoid, () )
    return
end

function call_genapwfr()
    ccall( (:my_genapwfr_, LIBLAPW), Cvoid, () )
    return
end

function call_genlofr()
    ccall( (:my_genlofr_, LIBLAPW), Cvoid, () )
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
    ccall( (:potks_, LIBLAPW), Cvoid, (Ref{Bool},), true )
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
