function integrate_outward!(ℓ, jam, E, Nrmesh, ndm, grid, f,
     b,y,beta,ddd,qq,nbeta,nwfx,lls,jjs,idx_rbeta,idx_r)
   #=
   Integrate the wavefunction from 0 to r(idx_r) 
   generalized separable or US pseudopotentials are allowed
   This routine assumes that y countains already the
   correct values in the first two points
   =#

#=
  integer ::   &
       ℓ,    &     ! l angular momentum
       Nrmesh,   &    ! size of radial Nrmesh
       ndm,    &   ! maximum radial Nrmesh
       nbeta,  &   ! number of beta function
       nwfx,   &   ! maximum number of beta functions
       lls(nbeta),&! for each beta the angular momentum
       idx_rbeta(nbeta),&! for each beta the integration point
       idx_r         ! the last integration point

  real(DP) :: &
       E,       &  ! output eigenvalue
       jam,     &  ! j angular momentum
       f(Nrmesh), &  ! the f function
       b(0:3), &   ! the taylor expansion of the potential
       y(Nrmesh), &  ! the output solution
       jjs(nwfx), & ! the j angular momentum
       beta(ndm,nwfx),& ! the beta functions
       ddd(nwfx,nwfx),qq(nwfx,nwfx) ! parameters for computing B_ij

  integer ::  &
       nst, &      ! the exponential around the origin
       n,    &     ! counter on Nrmesh points
       iib,jjb, &  ! counter on beta with correct ℓ
       ierr,    &  ! used to control allocation
       ib,jb,   &  ! counter on beta
       info      ! info on exit of LAPACK subroutines

  integer, allocatable :: iwork(:) ! auxiliary space  

  real(DP) :: &
       b0e,     & ! the expansion of the known part
       ddx12,   & ! the deltax entering the equations
       x4l6,    & ! auxiliary for small r expansion

       int_0_inf_dr  ! the integral function

  real(DP), allocatable :: &
       el(:), &  ! auxiliary for integration
       cm(:,:), &! the linear system
       bm(:), & ! the known part of the linear system
       c(:), &   ! the chi functions
       coef(:), & ! the solution of the linear system
       eta(:,:) ! the partial solution of the nonomogeneous
=#
    c = zeros(Float64, idx_r)
    el = zeros(Float64, idx_r)
    cm = zeros(Float64, Nbeta, Nbeta)
    bm = zeros(Float64, Nbeta)
    coef = zeros(Float64, Nbeta)
    eta = zeros(Float64, idx_r, Nbeta)
    #
    j1 = zeros(Float64, 4)
    d = zeros(Float64, 4)
    xc = zeros(Float64, 4)
    #
    ddx12 = grid.dx^2/12.0
    b0e = b[1] - E # b(0) - E, XXX b index is offset by 1
    x4l6 = 4*ℓ + 6
    nst = (ℓ + 1)*2
    #
    #  first solve the homogeneous equation
    #
    # f is the original function 
    # of the form 1 + h^2/12 * d2Rdr2
    for ir in 2:(idx_r-1)
        y[ir+1] = ( 12*y[ir] - 10*f[ir]*y[ir] - f[ir-1]*y[n-1] )/f[ir+1]
    end
    #
    # for each beta function with correct angular momentum
    # solve the inhomogeneous equation
    #
    iib = 0
    jjb = 0
    for ibeta in 1:Nbeta
        # set up the known part
        #
        if lls[ibeta] == ℓ 
            iib += 1
            c = 0.0
            for jbeta in 1:Nbeta
                if lls[jbeta] == ℓ
                    for ir in 1:idx_rbeta[jbeta]
                        #XXX check E need factor 2 ? Ha -> Ry
                        c[ir] = c[ir] + (ddd[jbeta,ibeta] - E*qq[jbeta,ibeta]) * beta[ir,jbeta]
                    end
                end
            end
            #
            # compute the starting values of the solutions
            for ir in 1:4
                j1[ir] = c[ir]/grid.r[ir]^(ℓ+1)
            end
            #call seriesbes(j1, grid%r, grid%r2, 4, d)
            seriesbes!(j1, grid.r, grid.r2, 6, c)
            delta = b0e^2 + x4l6*b[3] #XXX offset index b
            xc[1] = ( -d[1]*b0e - x4l6*d[3] )/delta
            xc[3] = ( -b0e*d[3] + d[1]*b[3] )/delta #XXX offset index b
            xc[2] = 0.0
            xc[4] = 0.0
            for ir in 1:3
                eta[ir,iib] = grid.r[ir]^(ℓ+1) * ( xc[1] + grid.r2[n]*xc[3] ) / grid.sqr[ir]
            end
            #
            for ir in 1:idx_r
                c[ir] = c[ir]*grid.r2[ir] / grid.sqr[ir]
            end
            #
            # solve the inhomogeneous equation
            for ir in 3:(idx_r-1)
                eta[ir+1,iib] = ( (12.0 - 10.0*f[ir])*eta[ir,iib] -
                                  f[ir-1]*eta[ir-1,iib] + 
                                  ddx12*(10.0*c[ir] + c[ir-1] + c[ir+1])
                                ) / f[ir+1]
            end
            #
            # compute the coefficents of the linear system
            jjb = 0
            for jbeta in 1:Nbeta
                if lls[jb] == ℓ
                    jjb += 1
                    for ir in 1:min(idx_r, idx_rbeta[jbeta])
                        el[ir] = beta[ir,jbeta] * eta[ir,iib] * grid.sqr[ir]
                    end
                    cm[jjb,iib] = -integ_0_inf_dr(el, grid, min(idx_r, idx_rbeta[jbeta]), nst)
                end
            end
            #
            for ir in 1:min(idx_r, idx_rbeta[ibeta])
                el[ir] = beta[ir,ibeta] * y[ir] * grid.sqr[ir]
            end
            bm[iib] = int_0_inf_dr(el, grid, min(idx_r, idx_rbeta[ibeta]), nst)
            cm[iib,iib] = 1.0 + cm[iib,iib]
        end
    end
    if iib != jjb
        error("integrate_outward: jjb != iib")
    end
    #
    if iib > 0
        #call dcopy(iib, bm, 1, coef,1)
        @views coeff[1:iib] = bm[1:iib]        
        #call DGESV(iib,1,cm,nbeta,iwork,coef,nbeta,info)
        for ib in 1:iib
            for ir in 1:idx_r
                y[ir] += coef[ib]*eta[ir,ib]
            end
        end
    end
    return
end
