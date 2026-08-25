function calc_pseudo_q!(ld1x_input, grid, qfunc, qfuncl, idx_rbeta)
    #=
  USE kinds, ONLY : DP
  USE ld1_parameters, ONLY : nwfsx
  USE ld1inc, ONLY : rcut, lls,  grid, ndmx, lmx2, nbeta, ikk, ecutrho, &
              rmatch_augfun, rmatch_augfun_nc
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: qfunc(ndmx,nwfsx,nwfsx)
  REAL(DP), INTENT(OUT) :: qfuncl(ndmx,nwfsx,nwfsx,0:lmx2)
  REAL(DP),  EXTERNAL :: int_0_inf_dr
  !
  ! variables for aug. functions generation
  ! 
  INTEGER  :: irc, ns, ns1, l1, l2, l3, lll, mesh, n, ik
  INTEGER  :: l1_e, l2_e
  REAL(DP) :: aux(ndmx)
  REAL(DP) :: augmom, ecutrhoq, rmatch
    =#
    
    Nrmesh = grid.Nrmesh
    Nbeta = ld1x_input.Nbeta
    lls = ld1x_input.lls
    rcut = ld1x_input.rcut # cutoff for projectors
    #
    ecutrho = 0.0
    l1_e = -1 # invalid value
    l2_e = -1 # invalid value
    fill!(qfuncl, 0.0)
    for ibeta in 1:Nbeta
        l1 = lls[ibeta]
        for jbeta in ibeta:Nbeta
            l2 = lls[jbeta]
            # Find the matching point
            idx_r = 0
            #if rmatch_augfun_nc
            # if `true` the norm conserving radii are used to pseudize the q functions
            rmatch = min(rcut[ibeta], rcut[jbeta])
            #else
            #rmatch = rmatch_augfun
            #end
            #
            for ir in 1:Nrmesh
                if grid.r[ir] > rmatch
                    idx_r = ir
                    break
                end
            end
            if (idx_r == 0) || (idx_r > Nrmesh-20)
                error("Something wrong with rmatch (too large?)")
            end
            #
            # Do the pseudization
            for l3 in range(abs(l1-l2), stop = (l1+l2), step = 2)
                @views ecutrhoq = compute_q_3bess!(grid, l3, l1+l2, idx_r, qfunc[:,ibeta,jbeta], qfuncl[:,ibeta,jbeta,l3])
                if ecutrhoq > ecutrho
                    ecutrho = ecutrhoq
                    l1_e = l1
                    l2_e = l2
                end
                qfuncl[1:Nrmesh,jbeta,ibeta,l3] = qfuncl[1:Nrmesh,ibeta,jbeta,l3]
            end
        end
    end
    #
    # Check that multipoles have not changed
    ir_c = maximum(idx_rbeta[1:Nbeta]) + 8
    augmom = 0.0
    aux = zeros(Float64, Nrmesh)
    for ibeta in 1:Nbeta
        l1 = lls[ibeta]
        for jbeta in ibeta:Nbeta
            l2 = lls[jbeta]
            for l3 in range(abs(l1-l2), stop = l1+l2, step = 2)
                for ir in 1:ir_c
                    aux[ir] = (qfuncl[ir,ibeta,jbeta,l3] - qfunc[ir,ibeta,jbeta]) * grid.r[ir]^l3
                end
                lll = l1 + l2 + 2 + l3
                augmom = integ_0_inf_dr(aux[1:ir_c], grid, ir_c,lll)
                if abs(augmom) > 1e-5
                    println("Problem with multipole ibeta=$ibeta,l1=$l1 jbeta=$jbeta,l2=$l2 l3=$l3 augmom=$augmom")
                end
            end
        end
    end
    println("Q pseudized with Bessel functions")
    println("Expected ecutrho= $ecutrho due to l1=$l1_e l2=$l2_e")

    return
end


# return ecutrho
function compute_q_3bess!(grid, ldip, ℓ, idx_r, chir, phi_out)
    #=
    This routine computes the phi_out function by pseudizing the
    chir function with a linear combination of three Bessel functions
    multiplied by r**2. In input it receives the point
    idx_r where the cut is done, the angular momentum lam of the 
    bessel functions and the function chir.
    Phi_out has the same ldip dipole moment of chir.
    =#

    #=

  integer ::    &
       ldip,    & ! input: the order of the dipole
       lam,     & ! input: the angular momentum
       idx_r         ! input: the point corresponding to rc

  real(DP) :: &
       xc(8)      ! output: the coefficients of the Bessel functions

  real(DP) ::         &
       chir(Nrmesh),    &   ! input: the all-electron function
       phi_out(Nrmesh)      ! output: the phi function
  !
  real(DP) ::  &
       ecutrho,& ! the expected cut-off on the charge density for this q
       fae,    & ! the value of the all-electron function
       f1ae,   & ! its first derivative
       f2ae,   & ! the second derivative
       dip       ! the norm of the function
    =#

    nbes = 3 #XXX HARCODED

    Nrmesh = grid.Nrmesh
    gi = zeros(Float64, Nrmesh)
    j1 = zeros(Float64, Nrmesh, nbes)
    cm = zeros(Float64, 3)
    bm = zeros(Float64, 3)
    xc = zeros(Float64, 8)

    nst = ℓ + 2 + ldip

    # compute the first and second derivative of input function at r(idx_r)
    fae = chir[idx_r]
    f1ae = deriv_7pts(chir, idx_r, grid.r[idx_r], grid.dx)
    f2ae = deriv2_7pts(chir, idx_r, grid.r[idx_r], grid.dx)
    #
    # compute the ldip dipole moment of the input function
    for ir in 1:(idx_r+1)
        gi[ir] = chir[ir] * grid.r[ir]^ldip  
    end
    dip = integ_0_inf_dr(gi, grid, idx_r, nst)
    #
    # RRKJ: the pseudo-wavefunction is written as an expansion into 3  
    #       spherical Bessel functions for r < r(idx_r)
    # find q_i with the correct log derivatives
    #call find_qi(f1ae/fae, xc[nbes+1:], idx_r, ldip, nbes, 2, iok)
    @views ld1x_find_qi!(grid, f1ae/fae, xc[(nbes+1):end], idx_r, ldip, nbes, 2) # flag is 2
    # iok is not used ? error is handled in ld1x_find_qi
    #
    # compute the Bessel functions and multiply by r**2
    for ibes in 1:nbes
        #call sph_bes(idx_r + 5, grid%r, xc(nbes+nc), ldip, j1(1,nc))
        for ir in 1:(idx_r+5)
            j1[ir,ibes] = sphericalbesselj(ldip, grid.r[ir]*xc[nbes+ibes])
        end
        jnor = j1[idx_r,ibes]*grid.r2[idx_r]
        for ir in 1:idx_r+5
            j1[ir,ibes] = j1[ir,ibes]*grid.r2[ir]*chir[idx_r]/jnor
        end
    end
    #
    # compute the bm functions (second derivative of the j1)
    # and the ldip dipole moment of the Bessel function (cm)
    for ibes in 1:nbes
        bm[ibes] = deriv2_7pts(j1[:,ibes], idx_r, grid.r[idx_r], grid.dx)
        for ir in 1:idx_r
            gi[ir] = j1[ir,ibes]*grid.r[ir]^ldip
        end
        cm[ibes] = integ_0_inf_dr(gi, grid, idx_r, nst)
    end
    #
    # solve the linear system to find the coefficients
    gam = ( bm[3] - bm[1] )/( bm[2] - bm[1] )
    delta = ( f2ae - bm[1] )/( bm[2] - bm[1] )
    #
    xc[3] = (dip - cm[1] + delta*(cm[1] - cm[2]))/( gam*(cm[1]-cm[2]) + cm[3] - cm[1] )
    xc[2] = -xc[3]*gam + delta
    xc[1] = 1.0 - xc[2] - xc[3]
    #
    # Set the function for r <= r[idx_r]
    for ir in 1:idx_r
        phi_out[ir] = xc[1]*j1[ir,1] + xc[2]*j1[ir,2] + xc[3]*j1[ir,3]
    end
    #
    # for r > r(idx_r) the function does not change
    for ir in (idx_r+1):Nrmesh
        phi_out[ir] = chir[ir]
    end
    ecutrho = 2.0*xc[6]^2
    return ecutrho
end