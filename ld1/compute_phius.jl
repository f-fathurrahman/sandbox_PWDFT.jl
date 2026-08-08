function compute_phius!(grid, ℓ, idx_rc, psi_in, phi_out, xc, iflag, els_in)
#=
This routine computes the phi functions by pseudizing the
all_electron chi functions. In input it receives, the point
ik where the cut is done, the angular momentum lam,
and the correspondence with the all eletron wavefunction
=#

    Nrmesh = grid.Nrmesh
    bm = zeros(Float64, 2)
    fact_nrm = zeros(Float64, 2)
    j1 = zeros(Float64, Nrmesh, 8)
    
    fill!(xc, 0.0) # reset xc
    #
    # compute first and second derivative
    #
    fae = psi_in[idx_rc]
    f1ae = deriv_7pts(psi_in, idx_rc, grid.r[idx_rc], grid.dx)
    f2ae = deriv2_7pts(psi_in, idx_rc, grid.r[idx_rc], grid.dx)
    #
    # find the q_i of the bessel functions
    #      
    #call find_qi(f1ae/fae,xc(4),ik,lam,2,1,iok)
    ncn = 2
    @views ld1x_find_qi!(grid, f1ae/fae, xc[4:4+ncn-1], idx_rc, ℓ, 2, 1) # or just use
    # Need to check iok?
    # compute the functions
    for ic in 1:2 # using two Bessel functions
        #call sph_bes(ik+5, grid%r, xc(3+nc), lam, j1(1,nc))
        for ir in 1:(idx_rc+5)
            j1[ir,ic] = sphericalbesselj(ℓ, xc[3+ic]*grid.r[ir])
        end
        fact_nrm[ic] = psi_in[idx_rc] / ( j1[idx_rc,ic] * grid.r[idx_rc] )
        for ir in 1:(idx_rc+5)
            j1[ir,ic] = j1[ir,ic] * grid.r[ir] * fact_nrm[ic]
        end
    end
    #
    # compute the second derivative and impose continuity of zero, 
    # first and second derivative
    for ic in 1:2
        bm[ic] = deriv2_7pts( j1[:,ic], idx_rc, grid.r[idx_rc], grid.dx )
    end
    xc[2] = (f2ae - bm[1])/(bm[2] - bm[1])
    xc[1] = 1.0 - xc[2]
    
    #if( iflag == 1 ) then
    #write(stdout,110) els_in,grid%r(ik),2.0_dp*xc(5)**2
    #110 format (5x, ' Wfc-us ',a3,' rcutus=',f6.3,'  Estimated cut-off energy= ', f8.2,' Ry')
    #if( verbosity == 'high') then
    #  write(stdout,'(5x,a)') 'rc*logder, xc(1), xc(2), rc*q(1),rc*q(2)'
    #  write(stdout,'(7f12.5)') grid%r(ik)*f1ae/fae, xc(1:2), xc(4:5)*grid%r(ik)
    #endif
    #endif
    #
    # define the phis function
    #
    for ir in 1:idx_rc
        phi_out[ir] = xc[1]*j1[ir,1] + xc[2]*j1[ir,2]
    end
    #
    for ir in (idx_rc+1):grid.Nrmesh
        phi_out[ir] = psi_in[ir]
    end
    #
    for ic in 1:2
        xc[ic] = xc[ic] * fact_nrm[ic]
    end
    #
    return
end
