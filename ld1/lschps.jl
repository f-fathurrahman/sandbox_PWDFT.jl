#=
Copyright (C) 2004-2010 Quantum ESPRESSO group
This file is distributed under the terms of the
GNU General Public License. See the file `License'
in the root directory of the present distribution,
or http://www.gnu.org/copyleft/gpl.txt .
=#

#
# adams extrapolation and interpolation formulas for
# outward and inward integration, abramowitz and stegun, p. 896
function lschps_aeo(y, j)
    # real(DP):: y(j), aeo
    res = 4.16666666667e-2 * ( 55.0*y[j] - 59.0*y[j-1] + 37.0*y[j-2] - 9.0*y[j-3] )
    return res
end


function lschps_aii(y, j)
    res = -4.16666666667e-2*( 9.0*y[j-1] + 19.0*y[j] - 5.0*y[j+1] + y[j+2] )
    return res
end

function lschps_eval_derv!(mmax, al, r, v, dv)
    # dv = dv/dr
    dv[1] = (-50*v[1] + 96*v[2] - 72*v[3] + 32*v(4) - 6*v[5] ) / ( 24*al*r[1] )
    dv[2] = (-6.0*v[1] - 20*v[2] + 36*v[3] - 12*v[4] + 2*v[5] ) / ( 24*al*r[2] )
    for i in 3:(mmax-2)
        dv[i] = (2*v[i-2] - 16*v[i-1] + 16*v[i+1] - 2*v[i+2] ) / ( 24*al*r[i] )
    end
    dv[mmax-1] = ( 3*v[mmax] + 10*v[mmax-1] - 18*v[mmax-2] + 6*v[mmax-3] - v[mmax-4] ) / ( 12*al*r[mmax-1] )
    dv[mmax] = ( 25*v[mmax] - 48*v[mmax-1] + 36*v[mmax-2] - 16*v[mmax-3] + 3*v[mmax-4] ) / ( 12*al*r[mmax] )
    return
end


#=
integrates radial pauli-type scalar-relativistic equation
on a logarithmic grid
modified routine to be used in finding norm-conserving
pseudopotential
  #
  # on input:
  #   mode = 1 find energy and wavefunction of bound states,
  #            scalar-relativistic (all-electron)
  #   mode = 2 find energy and wavefunction of bound state,
  #            nonrelativistic (pseudopotentials)
  #   mode = 3 fixed-energy calculation, for logarithmic derivatives
  #   mode = 4 find energy which produces a specified logarithmic
  #            derivative (nonrelativistic, pseudopotentials)
  #   mode = 5 is for pseudopotential to produce wavefunction beyond
  #            radius used for pseudopotential construction
  #   z    = atomic number
  #   eps  = convergence factor: eiganvalue is considered converged if
  #          the correction to eigenvalue is smaller in magnitude than
  #          eps times the magnitude of the current guess
  #   grid = structure containing radial grid information
  #   l, n = main and angular quantum numbers
  #   e    = starting estimate of the energy (mode=1,2)
  #          fixed energy at which the wavefctn is calculated (mode=3,4)
  #   v(i) = self-consistent potential
  #   nin  = integration up to r(nin) (mode=3,4,5)
  #
  # on output:
  #   e    = final energy (mode=1,2)
  #   u(i) = radial wavefunction (defined as the radial part of the wavefct
  #          multiplied by r)
  #   nstop= 0 if regular termination, 1 otherwise
  #   nin  = last grid point for which the wavefct is calculated (mode=1,2)
=#

function lschps( mode, z, eps, grid, nin, n, l, e, v, u, nstop )

#=
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: radial_grid_type
  USE ld1inc, ONLY : cau_fact
  IMPLICIT NONE
  #
  # I/O variables
  #
  INTEGER, INTENT (in) :: mode, n, l
  REAL(DP), INTENT(in) :: z, eps
  TYPE (radial_grid_type), INTENT(in) :: grid
  REAL(DP), INTENT(in) :: v(grid%mesh)
  INTEGER, INTENT(inout) :: nin
  REAL(DP), INTENT(inout) :: e
  INTEGER, INTENT(out) :: nstop
  REAL(DP), INTENT(out) :: u(grid%mesh)
  #
  # local variables
  #
  INTEGER, PARAMETER :: maxter=60
  REAL(DP), EXTERNAL:: aei, aeo, aii, aio
  # arrays  used as work space
  REAL(DP),ALLOCATABLE :: up(:),upp(:),cf(:),dv(:),fr(:),frp(:)
  REAL(DP):: al, als, cn
  real(DP):: de, emax, emin
  real(DP):: fss, gamma, ro, sc
  real(DP):: sls, sn, uld, uout,  upin, upout
  real(DP):: xkap
  INTEGER:: i, it, mmax, n_it, node, mch, ierr
=#
    nstop = 0
    al   = grid.dx
    mmax = grid.mesh

    up = zeros(Float64, mmax)
    upp = zeros(Float64, mmax)
    cf = zeros(Float64, mmax)
    dv = zeros(Float64, mmax)
    fr = zeros(Float64, mmax)
    frp = zeros(Float64, mmax)

    uld = 0.0
    if mode == 1 || mode == 3
        # relativistic calculation
        # fss=(1.0_dp/137.036_dp)^2
        fss = (1.0/cau_fact)^2
        if l == 0
           γ = sqrt(1.0 - fss*z^2)
        else
           γ = ( l*sqrt(l^2 - fss*z^2) + (l+1) * sqrt((l+1)^2 - fss*z^2) ) / (2*l + 1)
        end
    else
        # non-relativistic calculation
        fss = 1.0e-20
        γ = l + 1
    end
   
    sls = l*(l + 1)
    #
    # emin, emax = estimated bounds for e
    #
    if (mode == 1) || (mode == 2)
        emax = v[mmax] + sls/grid.r[mmax]^2
        emin = 0.0
        for i in 1:mmax
           emin = min(emin, v[i] + sls/grid.r[i]^2)
        end
        if e > emax
           e =1.25*emax
        end
        if e < emin
           e = 0.75*emin
        end
        if e > emax
           e = 0.5*(emax+emin)
        end
    elseif mode == 4
        emax = e + 10.0
        emin = e - 10.0
    end
    #
    for i in 1:4
        u[i] = 0.0
        up[i] = 0.0
        upp[i] = 0.0
    end
    als = al^2
    #
    # calculate dv/dr for darwin correction
    #
    lschps_eval_derv!(mmax, al, grid%r, v, dv )
    #
    # starting of loop on energy for bound state
    #
    for n_it in 1:maxter
        #
        # coefficient array for u in differential eq.
        for i in 1:mmax
           cf[i] = als*(sls + ( v[i] - e)*grid.r[i]^2)
        end
        #
        # find classical turning point for matching
        #
        if (mode == 1) || (mode == 2)
           for i in range(mmax, stop=2, step=-1)
              if ( cf[i-1] <= 0.0 ) && ( cf[i] > 0.0 )
                 mch = i
                 @goto LABEL40
              end
           end
           #PRINT '('' warning: wfc '',2i2,'' no turning point'')', n, l
           e = 0.0
           for i in 1:mmax
              u[i] = 0.0
           end
           nstop = 1
           @goto LABEL999
        else
           mch = nin
        end
      
        @label LABEL40
        # relativistic coefficient arrays for u (fr) and up (frp).
        for i in 1:mmax
            fr[i] = als*( grid.r[i]^2 )*0.25 * (
               -fss*(v[i] - e)^2 + fss*dv[i] / ( grid.r[i] * ( 1.0 + 0.25*fss*(e - v[i]) ) )
            )
            frp[i] = -al*grid.r[i]*0.25 * fss*dv[i] / ( 1.0 + 0.25 * fss*(e - v[i]) )
        end
        #
        # start wavefunction with series
        #
        for i in 1,4
            u[i] = grid.r[i]^γ
            up[i] = al*γ*grid.r[i]^γ
            upp[i] = (al + frp[i]) * up[i] + ( cf[i] + fr[i] ) * u[i]
        end
        #
        # outward integration using predictor once, corrector
        # twice
        node = 0
        #
        for i in 4:(mch-1)
            u[i+1] = u[i] + lschps_aeo(up, i)
            up[i+1] = up[i] + lschps_aeo(upp, i)
            for it in 1:2
                upp[i+1] = (al + frp[i+1] )*up[i+1] + ( cf[i+1] + fr[i+1] ) * u[i+1]
                up[i+1] = up[i] + lschps_aio(upp, i)
                u[i+1] = u[i] + lschps_aio(up, i)
            end
            if u[i+1]*u[i] <= 0.0
                node = node + 1
            end
            #
            uout = u[mch]
            upout = up[mch]
            #
            if (node-n+l+1 == 0) || (mode == 3) || (mode == 5)
            #
            if (mode == 1) || (mode == 2)
                #
                # start inward integration at 10*classical turning
                # point with simple exponential
                nin = mch + 2.3/al
                if (nin + 4) > mmax
                    nin = mmax - 4
                end
                xkap = sqrt( sls/grid.r[nin]^2 + 2.0*(v[nin] - e) )
                #
                for i in nin:(nin+4)
                    u[i] = exp( -xkap * ( grid.r[i] - grid.r[nin] ) )
                    up[i] = -grid.r[i] * al * xkap * u[i]
                    upp[i] = (al + frp[i] ) * up[i] + (cf[i] + fr[i] ) * u[i]
                end
                #
                # integrate inward
                #
                for i in range(nin, stop=mch+1, step=-1)
                    u[i-1] = u[i] + lschps_aei(up, i)
                    up[i-1] = up[i] + lschps_aei(upp, i)
                    for it in 1:2
                        upp[i-1] = (al + frp[i-1] )*up[i-1] + (cf[i-1] + fr[i-1] )*u[i-1]
                        up[i-1] = up[i] + lschps_aii(upp, i)
                        u[i-1] = u[i] + lschps_aii(up, i)
                    end
                end
                #
                # scale outside wf for continuity
                sc = uout/u[mch]
                #
                for i in mch:nin
                    up[i] = sc*up[i]
                    u[i] = sc*u[i]
                end
                #
                upin = up[mch]
            #
            else
                #
                upin = uld*uout
                #
            end # if
            #
            # perform normalization sum
            #
            ro = grid.r[1] * exp(-0.5 * grid.dx)
            sn = ro^(2*γ + 1)/(2*γ + 1)
            #
            for i in 1:(nin-3)
                sn = sn + al*grid.r[i] * u[i]^2
            end
            #
            sn = sn + al*(23.0*grid.r[nin-2]*u[nin-2]^2 + 28.0*grid.r[nin-1] * u[nin-1]^2 +
                          9.0*grid.r[nin] * u[nin]^2) / 24.0
            #
            # normalize u
            cn = 1.0/sqrt(sn)
            uout = cn*uout
            upout = cn*upout
            upin = cn*upin
            #
            for i in 1:nin
                up[i] = cn*up[i]
                u[i] = cn*u[i]
            end
            for i in (nin+1):mmax
                u[i] = 0.0
            end
            #
            # exit for fixed-energy calculation
            #
            if (mode == 3) || (mode == 5)
                @goto LABEL999
            end
            # perturbation theory for energy shift
            de = uout*(upout - upin)/( al*grid.r[mch] )
            #
            # convergence test and possible exit
            #
            if abs(de) < max(abs(e), 0.2)*eps
                @goto LABEL999
            end
            #
            if de > 0.0
                emin = e
            else
                emax = e
            end
            e = e + de
            if (e > emax) || (e < emin)
                e = 0.5*( emax + emin )
            end
            #
        elseif node-n+l+1 < 0
            # too few nodes
            emin = e
            e = 0.5*(emin + emax)
        else
            # too many nodes
            emax = e
            e = 0.5*(emin + emax)
        end # if
    end # do

    #PRINT '('' warning: wfc '',2i2,'' not converged'')', n, l
    u = 0.0
    nstop = 1

    @label LABEL999

    return

end