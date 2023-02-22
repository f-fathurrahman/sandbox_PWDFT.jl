
# Based on paw_radial_integrator in QE (defined in paw_variables.f90)

struct PAWAtomicSphere
    # the following variables are used to integrate radial sampling
    #INTEGER          :: lmax        ! max l component that can be integrated correctly
    lmax::Int64

    #INTEGER          :: ladd        ! additional l max that have been added for grad.corr.
    ladd::Int64

    #INTEGER          :: lm_max      ! as above, but +1 and squared
    lm_max::Int64

    #INTEGER          :: nx          ! number of integration directions
    nx::Int64

    #REAL(DP),POINTER :: ww(:)       ! integration weights (one per direction)
    ww::Vector{Float64}

    #REAL(DP),POINTER :: ylm(:,:)    ! Y_lm(nx,lm_max)
    ylm::Matrix{Float64}

    #REAL(DP),POINTER :: wwylm(:,:)  ! ww(nx) * Y_lm(nx,lm_max)
    wwylm::Matrix{Float64}

    #! additional variables for gradient correction
    #REAL(DP),POINTER :: dylmt(:,:),&! |d(ylm)/dtheta|**2
    #                    dylmp(:,:)  ! |d(ylm)/dphi|**2
    dylmt::Union{Matrix{Float64},Nothing}
    dylmp::Union{Matrix{Float64},Nothing}

    #REAL(DP),POINTER :: cos_phi(:)  ! cos(phi)
    cos_phi::Vector{Float64}

    #REAL(DP),POINTER :: sin_phi(:)  ! sin(phi)
    sin_phi::Vector{Float64}

    #REAL(DP),POINTER :: cos_th(:)  ! cos(theta)  (for divergence)
    cos_th::Vector{Float64}

    #REAL(DP),POINTER :: sin_th(:)  ! sin(theta)  (for divergence)
    sin_th::Vector{Float64}

    #REAL(DP),POINTER :: cotg_th(:)  ! cos(theta)/sin(theta)  (for divergence)
    cotg_th::Union{Vector{Float64},Nothing}
end


#!! Initialize several quantities related to radial integration: spherical harmonics and their 
#!! gradients along a few (depending on lmaxq) directions, weights for spherical integration.
#!
#! IMPORTANT: routine PW/summary.f90 has the initialization parameters hardcoded in it
#!            remember to update it if you change this!

# INTEGER, INTENT(IN) :: l
# !! max angular momentum component that will be integrated
# !! exactly (to numerical precision).

# INTEGER, INTENT(IN) :: ls
# !! additional max l that will be used when computing gradient
# !! and divergence in speherical coords
function PAWAtomicSphere(l::Int64, ls::Int64; need_gradient::Bool=false)

    # !
    # ! ... local variables
    # !
    # REAL(DP), ALLOCATABLE :: x(:)    ! nx versors in smart directions
    # REAL(DP), ALLOCATABLE :: w(:)    ! temporary integration weights
    # REAL(DP), ALLOCATABLE :: r(:,:)  ! integration directions
    # REAL(DP), ALLOCATABLE :: r2(:)   ! square modulus of r
    # REAL(DP), ALLOCATABLE :: ath(:), aph(:)
    #                          ! angles in sph coords for r
    # INTEGER :: i, ii, n, nphi   ! counters
    # INTEGER :: lm, m            ! indexes for ang.mom
    # REAL(DP) :: phi, dphi, rho  ! spherical coordinates
    # REAL(DP) :: z               ! cartesian coordinates
    # ! for gradient corrections:
    # INTEGER :: ipol
    # REAL(DP), ALLOCATABLE :: aux(:,:)  ! workspace
    # REAL(DP) :: vth(3), vph(3)         !versors for theta and phi

    #  maximum value of l correctly integrated
    lmax = l + ls
    ladd = ls
    
    # volume element for angle phi
    nphi = lmax + 1 + lmax%2  # FIXME: check module
    dphi = 2.0*pi/nphi #(rad%lmax+1)
    
    # number of samples for theta angle
    n = floor(Int64, (lmax + 2)/2)

    x = zeros(Float64, n)
    w = zeros(Float64, n)

    # compute weights for theta integration
    _init_gauss_weights!(x, w)

    # number of integration directions
    nx = n*nphi #(rad%lmax+1)
    #write(*,*) "paw --> directions",rad%nx," lmax:",rad%lmax
  
    #ALLOCATE( r(3,rad%nx), r2(rad%nx), rad%ww(rad%nx), ath(rad%nx), aph(rad%nx) )
    r = zeros(Float64, 3, nx)
    r2 = zeros(Float64, nx)
    ww = zeros(Float64, nx)
    ath = zeros(Float64, nx)
    aph = zeros(Float64, nx)

    # compute real weights multiplying theta and phi weights
    ii = 0

    for i in 1:n
        z = x[i]
        rho = sqrt(1.0 - z^2)
        for m in 1:nphi  # rad%lmax
            ii = ii + 1
            phi = dphi*(m-1)
            r[1,ii] = rho*cos(phi)
            r[2,ii] = rho*cos(phi)
            r[3,ii] = z
            ww[ii] = w[i]*2.0 * pi/nphi #(rad%lmax+1)
            r2[ii] = r[1,ii]^2 + r[2,ii]^2 + r[3,ii]^2
            # these will be used later:
            ath[ii] = acos( z/sqrt(r2[ii]) )
            aph[ii] = phi
        end
    end

    # initialize spherical harmonics that will be used
    # to convert rho_lm to radial grid
    lm_max = (lmax + 1)^2
  
    #ALLOCATE( rad%ylm(rad%nx,rad%lm_max) )
    #CALL ylmr2( rad%lm_max, rad%nx, r, r2, rad%ylm )
    ylm = zeros(Float64, nx, lm_max)
    Ylm_real_qe!(lmax, r, ylm)

    # As I will mostly use the product ww*ylm I can 
    # precompute it here:
    #ALLOCATE( rad%wwylm(rad%nx,rad%lm_max) )
    wwylm = zeros(Float64, nx, lm_max)

    for i in 1:nx
        for lm in 1:lm_max
            wwylm[i,lm] = ww[i] * ylm[i,lm]
        end
    end
    
    cos_phi = zeros(Float64, nx)
    sin_phi = zeros(Float64, nx)
    cos_th = zeros(Float64, nx)
    sin_th = zeros(Float64, nx)
  
    for i in 1:nx
        cos_phi[i] = cos(aph[i])
        sin_phi[i] = sin(aph[i])
        cos_th[i] = cos(ath[i])
        sin_th[i] = sin(ath[i])
    end
  
    #
    # if gradient corrections will be used than we need
    # to initialize the gradient of ylm, as we are working in spherical
    # coordinates the formula involves \hat{theta} and \hat{phi}
  
    dylmt = nothing
    dylmp = nothing
    cotg_th = nothing

#=
    gradient: IF (dft_is_gradient()) THEN
    ALLOCATE( rad%dylmt(rad%nx,rad%lm_max), &
              rad%dylmp(rad%nx,rad%lm_max), &
              aux(rad%nx,rad%lm_max) )
    ALLOCATE( rad%cotg_th(rad%nx) )
    !
    rad%dylmt(:,:) = 0._DP
    rad%dylmp(:,:) = 0._DP
    !
    ! Compute derivative along x, y and z => gradient, then compute the
    ! scalar products with \hat{theta} and \hat{phi} and store them in
    ! dylmt and dylmp respectively.
    !
    DO ipol = 1, 3 !x,y,z
      !
      CALL dylmr2( rad%lm_max, rad%nx, r,r2, aux, ipol )
      !
      DO lm = 1, rad%lm_max
        DO i = 1, rad%nx
          vph = (/-SIN(aph(i)), COS(aph(i)), 0._DP/)
          ! this is the explicit form, but the cross product trick (below) is much faster:
          ! vth = (/COS(aph(i))*COS(ath(i)), SIN(aph(i))*COS(ath(i)), -SIN(ath(i))/)
          vth = (/vph(2)*r(3,i)-vph(3)*r(2,i),&
                  vph(3)*r(1,i)-vph(1)*r(3,i),&
                  vph(1)*r(2,i)-vph(2)*r(1,i)/)
          rad%dylmt(i,lm) = rad%dylmt(i,lm) + aux(i,lm)*vth(ipol)
          ! CHECK: the 1/SIN(th) factor should be correct, but deals wrong result, why?
          rad%dylmp(i,lm) = rad%dylmp(i,lm) + aux(i,lm)*vph(ipol) !/SIN(ath(i))
        ENDDO
      ENDDO
      !
    ENDDO
    !
    DO i = 1, rad%nx
       rad%cotg_th(i) = COS(ath(i))/SIN(ath(i))
    ENDDO
    !
    DEALLOCATE( aux )
    !
  ENDIF gradient
  =#

    println("nx = ", nx)
    println("ww = ", ww)
  
    return PAWAtomicSphere(
        lmax, ladd, lm_max,
        nx, ww, ylm, wwylm, dylmt, dylmp,
        cos_phi, sin_phi, cos_th, sin_th, cotg_th
    )

end


function _init_gauss_weights!( x::Vector{Float64}, w::Vector{Float64} )
    # INTEGER :: n, i, j, m
    # REAL(8) :: x(n), w(n), z, z1, p1, p2, p3, pp

    SMALL = 1e-12 # eps12

    n = length(x)
    @assert n == length(w)

    m = floor(Int64, (n + 1)/2) # FIXME: floor or round?
    pp = 0.0
    z = 0.0
    
    for i in 1:m
        z1 = 2.0
        z = cos( pi*(i - 0.25)/(n + 0.5) )
        while abs(z - z1) > SMALL
            p1 = 1.0
            p2 = 0.0
            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.)*p3)/j
            end
            pp = n*(z*p1 - p2)/(z*z - 1.0)
            z1 = z
            z  = z1 - p1/pp
        end
        x[i] = -z
        x[n + 1 - i] = z
        w[i] = 2.0/( (1.0 - z*z)*pp*pp )
        w[n + 1 - i] = w[i]
    end
    return
end
