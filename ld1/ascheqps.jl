function ascheqps!( nam, ℓ, jam, E0, Nrmesh, ndm, grid, Vpot, thresh,
                    y, beta, ddd, qq, nbeta, nwfx, lls, jjs, ikk, nstop )
#=
numerical integration of a generalized radial schroedinger equation,
using Numerov with outward and inward integration and matching.
Works for both norm-conserving nonlocal and US pseudopotentials
Requires in input a good estimate "E0" of the energy

  integer, intent(in) :: &
       nam, &
       ℓ, &      ! l angular momentum
       Nrmesh,&      ! size of radial Nrmesh
       ndm, &      ! maximum radial Nrmesh 
       nbeta,&     ! number of beta function  
       nwfx, &     ! maximum number of beta functions
       ikk(nbeta),&! for each beta the point where it become zero
       lls(nbeta)  ! for each beta the angular momentum

  real(DP), intent(in) :: &
       jam,       & ! j angular momentum
       Vpot(Nrmesh),& ! the local potential 
       thresh,    & ! precision of eigenvalue
       jjs(nwfx), & ! the j angular momentum
       beta(ndm,nwfx), &            ! the beta functions
       ddd(nwfx,nwfx),qq(nwfx,nwfx) ! parameters for computing B_ij

  real(DP), intent(inout) :: &
       E0,      &  ! output eigenvalue
       y(Nrmesh)     ! the output solution

  integer, intent(out) :: &
       nstop       ! error code, used to check the behavior of the routine
  !
  !    the local variables
  !
  integer :: &
       ndcr,  &    ! number of required nodes
       n1, n2, &   ! counters
       ikl         ! auxiliary variables
  real(DP) :: &
       work(nbeta),& ! auxiliary space
       E,          &  ! energy
       ddx12,      &  ! dx^2/12 used for Numerov integration
       sqlhf,      &  ! the term for angular momentum in equation
       ze2,        &  ! possible coulomb term aroun the origin (set 0)
       b(0:3),     &  ! coefficients of taylor expansion of potential
       eup,elw,    & ! actual energy interval
       ymx,        & ! the maximum value of the function
       fe,integ,dfe,de, &! auxiliary for numerov computation of E
       eps,        & ! the epsilon of the delta E
       yln, xp, expn,& ! used to compute the tail of the solution
       int_0_inf_dr  ! integral function

  real(DP), allocatable :: &
       fun(:),  &   ! integrand function
       f(:),    &   ! the f function
       el(:),c(:) ! auxiliary for inward integration

  integer, parameter :: &
       NmaxIter=100    ! maximum number of iterations

  integer :: &
       n,  &    ! counter on Nrmesh points
       iterSch,&   ! counter on iteration
       idx_r,  &   ! matching point
       ns,  &   ! counter on beta functions
       l1,  &   ! ℓ+1
       nst, &   ! used in the integration routine
       ierr, &
       ncross,& ! actual number of nodes
       ir_start  ! starting point for inward integration

  logical, save :: first(0:10,0:10) = .true.
=#

    Nrmesh = grid.Nrmesh

    # set up constants and allocate variables the 
    fun = zeros(Float64, Nrmesh)
    f = zeros(Float64, Nrmesh)
    el = zeros(Float64, Nrmesh)
    c = zeros(Float64, Nrmesh)
    
    nstop = 0
    ir_start = 0
    E = E0
    # write(6,*) 'entering ', nam,ℓ, E
    eup = 0.3*E
    elw = 1.3*E
    ndcr = nam - ℓ - 1
    # println("entering ascheqps ", Vpot(Nrmesh-20)*grid%r(Nrmesh-20))

    ddx12 = grid.dx^2/12.0
    l1 = ℓ + 1
    nst = l1*2
    sqlhf = (ℓ + 0.5)^2
    #
    # series developement of the potential near the origin
    for ir in 1:4
        y[ir] = Vpot[ir]
    end
    radial_grid_series!( y, grid.r, grid.r2, b )
    println("enter ℓ=$ℓ, eup=$eup, elw=$elw, E=$E")
    #
    #  set up the f-function and determine the position of its last
    #  change of sign
    #  f < 0 (approximatively) means classically allowed   region
    #  f > 0         "           "        "      forbidden   "
    #
    for iterSch in 1:NmaxIter
        println("starting iterSch=$iterSch, elw=$elw, E=$E, eup=$eup")
        idx_r = 1
        f[1] = ddx12*(grid.r2[1]*(Vpot[1] - E) + sqlhf) # XXX change to Ha
        for ir in 2:Nrmesh
            f[ir] = ddx12*(grid.r2[ir]*(Vpot[ir] - E) + sqlhf)
            if ( f[i] != abs(f[i])*sign(f[i-1]) ) && (ir < Nrmesh-5)
                idx_r = ir
            end
        end
        if (idx_r == 1) || grid.r[idx_r] > 4.0
            idx_r = round(Int64, Nrmesh*3/4)
        end
        #
        if idx_r >= Nrmesh-2
            error("No point found for matching")
            # probably something wrong with the potential
        end
        #
        # determine if idx_r is sufficiently large
        for ibeta in 1:Nbeta
            if (lls[ibeta] == ℓ) && ( idx_rbeta[ibeta] > idx_r )
                idx_r = idx_rbeta[ibeta] + 3
            end
        end
        #
        # if everything is ok continue the integration and define f
        for ir in 1:Nrmesh
            f[ir] = 1.0 - f[ir]
        end
        #
        # determination of the wave-function in the first two points by
        # series developement
        #
        # no coulomb divergence in the origin for a pseudopotential
        ze2 = 0.0 
        start_scheq!( ℓ, E, b, grid, ze2, y )
        #
        # outward integration before idx_r
        #
        integrate_outward( ℓ, jam, E, Nrmesh, ndm, grid, f, b, y, beta, ddd, qq,
                          nbeta, nwfx, lls, jjs, idx_rbeta, idx_r)

        ncross = 0
        ymx = 0.0
        for ir in 2:(idx_r-1)
            # XXX why not simply compare the sign?
            if y[ir] != abs(y[ir])*sign(y[ir+1]) # XXX this is float comparison
                ncross += 1
            end
            ymx = max(ymx, abs(y[ir+1]))
        end
        # If at this point the number of nodes is wrong it means that something
        # is probably wrong in the calling routines. A ghost might be present
        # in the pseudopotential. With a nonlocal pseudopotential there is no
        # node theorem so strictly speaking the following instructions are
        # wrong but sometimes they help so we keep them here.
        if ndcr != ncross # `first` should be inout argument
            println("Warning: n=$nam l=$ℓ, expecting $ndcr nodes, found $ncross")
            println("Setting wfc to zero for this iteration")
        end
        #
        if ndcr < ncross
            # too many crossings. E is an upper bound to the true eigenvalue.
            # increase abs(E)
            #
            eup = E
            E = 0.9*elw + 0.1*eup
            println("too many crossing: ncross=$ncross, ndcr=$ndcr")
            fill!(y, 0.0)
            ymx = 0.0
            @goto LABEL300
        elseif ndcr > ncross
            #
            # too few crossings. E is a lower bound to the true eigenvalue.
            # decrease abs(E)
            #
            elw = E
            E = 0.9*eup + 0.1*elw
            println("Too few crossing ncross=$ncross, ndcr=$ndcr")
            fill!(y, 0.0) #XXX Need this? 
            ymx = 0.0
            @goto LABEL300
        end
        #
        # inward integration up to idx_r
        #
        integrate_inward!(E, Nrmesh, ndm, grid, f, y, c, el, idx_r, ir_start)
        #
        # if necessary, improve the trial eigenvalue by the cooley's procedure.
        # jw cooley math of comp 15,363(1961)
        #
        fe = (12.0 - 10.0*f[idx_r])*y[idx_r] - f[idx_r-1]*y[idx_r-1] - f[idx_r+1]*y[idx_r+1]
        #
        # adjust the normalization if needed
        if ymx >= 1.0e10
            @views y[:] = y[:]/ymx
        end
        #
        #  calculate the normalization
        for ibeta in 1:Nbeta
            if (ℓ == lls[ibeta]) # also need to check jj for relativistic case
                idx_r_l = idx_rbeta[ibeta]
                for ir in 1:idx_r_l
                    fun[ir] = beta[ir,jbeta]*y[ir]*grid.sqr[ir]
                end
                work[ibeta] = integ_0_inf_dr(fun, grid, ikl, nst)
            else
                work[ibeta] = 0.0
            end
        end
        #
        for ir in 1:ir_start
            fun[ir] = y[ir]*y[ir]*grid.r[ir]
        end
        ss = integ_0_inf_dr(fun, grid, ir_start, nst)
        for ibeta in 1:Nbeta, jbeta in 1:Nbeta
            ss += qq[ibeta,jbeta]*work[ibeta]*work[jbeta]
        end
        dfe = -y[idx_r]*f[idx_r]/grid.dx/integ
        de = -fe*dfe
        epsE = abs(de/E)
        #  write(6,'(i5, 3f20.12)') iterSch, E, de
        if abs(de) < thresh
            @goto LABEL600
        end
        #
        if epsE > 0.25
            de = 0.25*de/epsE
        end
        #
        if de > 0.0
            elw = E
        end
        if de < 0.0
            eup = E
        end
        E = E + de
        #
        if E > eup
            E = 0.9*eup + 0.1*elw
        end
        #
        if E < elw
            E = 0.9*elw + 0.1*eup
        end
        @LABEL300 continue
    end
    nstop = 1
    
    if ir_start == 0
        @goto LABEL900 # return?
    end
  
    @LABEL600
    #  
    # exponential tail of the solution if it was not computed
    #
    if ir_start < Nrmesh
        for ir in ir_start:Nrmesh-1
            if y[ir] == 0.0 # XXX Float comparison
                y[ir+1] = 0.0
            else
                yln = log(abs(y[ir]))
                xp = -sqrt(12.0*abs(1.0-f[ir]))
                expn = yln + xp
                if expn < -80.0
                    y[ir+1] = 0.0
                else
                    y[ir+1] = abs(exp(expn))*sign(y[ir])
                end
            end
        end
    end
    #
    # normalize the eigenfunction as if they were norm conserving. 
    # If this is a US PP the correct normalization is done outside this routine.
    #
    for ir in 1:Nrmesh
        el[ir] = grid.r[ir]*y[ir]*y[ir]
    end
    ss = int_0_inf_dr(el, grid, Nrmesh, nst)
    if ss > 0.0
        ss = sqrt(ss)
        for ir in 1:Nrmesh
            y[ir] = grid.sqr[ir] * y[ir] / ss
        end
        E0 = E
    else
        nstop = 1
    end
    @LABEL900
    return

end


#=

!--------------------------------------------------------------------------
subroutine my_ascheqps_drv(veff, ncom, thresh, flag_all, nerr)
!--------------------------------------------------------------------------

  ! This routine is a driver that calculates for the test
  ! configuration the solutions of the Kohn and Sham equation
  ! with a fixed pseudo-potential. The potentials are assumed
  ! to be screened. The effective potential veff is given in input.
  ! The output wavefunctions are written in phits and are normalized.
  ! If flag is .true. compute all wavefunctions, otherwise only
  ! the wavefunctions with positive occupation.
  !      
  use kinds, only: dp
  use ld1_parameters, only: nwfsx
  use radial_grids, only: ndmx
  use ld1inc, only: grid, pseudotype, rel, &
                    lls, jjs, qq, ikk, ddd, betas, nbeta, vnl, &
                    nwfts, iswts, octs, llts, jjts, nnts, enlts, phits 
  implicit none

  integer ::    &
          nerr, &     ! control the errors of the routine ascheqps
          ncom        ! number of components of the pseudopotential

  real(DP) :: &
       veff(ndmx,ncom)    ! work space for writing the potential 

  logical :: flag_all    ! if true calculates all the wavefunctions

  integer ::  &
       ns,    &  ! counter on pseudo functions
       is,    &  ! counter on spin
       nbf,   &  ! auxiliary nbeta
       n,     &  ! index on r point
       nstop, &  ! errors in each wavefunction
       ind

  real(DP) :: &
       vaux(ndmx,2)     ! work space for writing the potential 

  real(DP) :: thresh         ! threshold for selfconsistency
  
  write(*,*)
  write(*,*) '<div> ENTER my_ascheqps_drv'
  write(*,*)
  
  !
  ! compute the pseudowavefunctions in the test configuration
  !
  if (pseudotype == 1) then
    nbf = 0
  else
    nbf = nbeta
  endif

  nerr = 0
  do ns = 1,nwfts
    if( octs(ns) > 0.0 .or. ( octs(ns) > -1.0 .and. flag_all ) ) then
      is = iswts(ns)
      if( ncom==1 .and. is==2) then
        call errore('ascheqps_drv','incompatible spin',1)
      endif
      !
      if( pseudotype == 1 ) then
        !
        if( rel < 2 .or. llts(ns) == 0 .or. &
          & abs(jjts(ns)-llts(ns)+0.5) < 0.001) then
          ind = 1
        !
        elseif( rel == 2 .and. llts(ns) > 0 .and. &
              & abs(jjts(ns)-llts(ns)-0.5) < 0.001) then
          ind = 2
        else
          call errore('my_ascheqps_drv', 'unexpected case', 1)
        endif
        !
        do n = 1,grid%Nrmesh
          vaux(n,is) = veff(n,is) + vnl(n,llts(ns),ind)
        enddo
      else
        ! other pseudotypes
        do n = 1,grid%Nrmesh
          vaux(n,is) = veff(n,is)
        enddo
      endif
      !
      call my_ascheqps( nnts(ns),llts(ns),jjts(ns),enlts(ns),grid%Nrmesh,ndmx,&
                    &   grid,vaux(1,is),thresh,phits(1,ns),betas,ddd(1,1,is),qq,nbf, &
                    &   nwfsx,lls,jjs,ikk,nstop)
      write(*,*) ns, nnts(ns),llts(ns), jjts(ns), enlts(ns)
      !
      ! normalize the wavefunctions 
      !
      call normalize(phits(1,ns), llts(ns), jjts(ns), ns)
      !
      !   not sure whether the "best" error code should be like this:
      ! IF ( octs(ns) > 0.0 ) nerr = nerr + nstop
      !   i.E. only for occupied states, or like this:
      nerr = nerr + nstop
    endif ! if octs is larger than zero
  enddo

  write(*,*)
  write(*,*) '</div> EXIT my_ascheqps_drv'
  write(*,*)

  return
end subroutine
=#

