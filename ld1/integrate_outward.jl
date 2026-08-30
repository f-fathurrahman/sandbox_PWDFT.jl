function integrate_outward!(ℓ, E, grid, f, b, y, beta, ddd, qq, lls, idx_rbeta, idx_r)
   #=
   Integrate the wavefunction from 0 to r(idx_r) 
   generalized separable or US pseudopotentials are allowed
   This routine assumes that y countains already the
   correct values in the first two points
   =#
    Nbeta = size(beta, 2)
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
    #b0e = 2*(b[1] - E) # b(0) - E, XXX b index is offset by 1
    b0e = b[1] - E # b(0) - E, XXX b index is offset by 1
    x4l6 = 4*ℓ + 6
    nst = (ℓ + 1)*2
    #
    #  first solve the homogeneous equation
    #
    # f is the original function 
    # of the form 1 + h^2/12 * d2Rdr2
    for ir in 2:(idx_r-1)
        y[ir+1] = ( 12*y[ir] - 10*f[ir]*y[ir] - f[ir-1]*y[ir-1] )/f[ir+1]
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
            fill!(c, 0.0)
            for jbeta in 1:Nbeta
                if lls[jbeta] == ℓ
                    for ir in 1:idx_rbeta[jbeta]
                        #XXX check E need factor 2 ? Ha -> Ry
                        c[ir] += ( ddd[jbeta,ibeta] - E*qq[jbeta,ibeta] ) * beta[ir,jbeta]
                        #c[ir] += 2 * ( ddd[jbeta,ibeta] - E*qq[jbeta,ibeta] ) * beta[ir,jbeta]
                    end
                end
            end
            #
            # compute the starting values of the solutions
            for ir in 1:4
                j1[ir] = c[ir]/grid.r[ir]^(ℓ+1)
            end
            #call seriesbes(j1, grid%r, grid%r2, 4, d)
            seriesbes!(j1, grid.r, grid.r2, 4, c)
            delta = b0e^2 + x4l6*b[3] #XXX offset index b
            xc[1] = ( -d[1]*b0e - x4l6*d[3] )/delta
            xc[3] = ( -b0e*d[3] + d[1]*b[3] )/delta #XXX offset index b
            xc[2] = 0.0
            xc[4] = 0.0
            for ir in 1:3
                eta[ir,iib] = grid.r[ir]^(ℓ+1) * ( xc[1] + grid.r2[ir]*xc[3] ) / sqrt(grid.r[ir])
            end
            #
            for ir in 1:idx_r
                c[ir] = c[ir]*grid.r2[ir] / sqrt(grid.r[ir])
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
                if lls[jbeta] == ℓ
                    jjb += 1
                    for ir in 1:min(idx_r, idx_rbeta[jbeta])
                        el[ir] = beta[ir,jbeta] * eta[ir,iib] * sqrt(grid.r[ir])
                    end
                    cm[jjb,iib] = -integ_0_inf_dr(el, grid, min(idx_r, idx_rbeta[jbeta]), nst)
                end
            end
            #
            for ir in 1:min(idx_r, idx_rbeta[ibeta])
                el[ir] = beta[ir,ibeta] * y[ir] * sqrt(grid.r[ir])
            end
            bm[iib] = integ_0_inf_dr(el, grid, min(idx_r, idx_rbeta[ibeta]), nst)
            cm[iib,iib] = 1.0 + cm[iib,iib]
        end
    end
    if iib != jjb
        error("integrate_outward: jjb != iib")
    end
    #
    if iib > 0
        #call dcopy(iib, bm, 1, coef,1)
        @views coef[1:iib] = bm[1:iib]        
        #call DGESV(iib,1,cm,nbeta,iwork,coef,nbeta,info)
        coef[1:iib] .= cm[1:iib,1:iib] \ coef[1:iib]
        for ib in 1:iib
            for ir in 1:idx_r
                y[ir] += coef[ib]*eta[ir,ib]
            end
        end
    end
    return
end
