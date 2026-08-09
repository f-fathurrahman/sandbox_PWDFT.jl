function _set_psi_in!(Zval, Vpot, grid, idx_rc, ℓ, E_in, psi_out)
    #=
    This subroutine calculates the all electron wavefunction psi at the
    input energy e
    =#

    #=
  USE kinds, ONLY : dp
  USE radial_grids, ONLY : ndmx
  USE ld1inc, ONLY : grid, rel, zed, vpot
  IMPLICIT NONE
  INTEGER :: l, ik    ! input: angular momentum and index of the cut-off radius
  REAL(DP) :: e, j    ! input: energy and total angular momentum
  REAL(DP) :: psi_out(ndmx) ! output: the function psi.
  REAL(DP) :: psi_out_rel(ndmx) ! output: the function psi (small component).
  REAL(DP) :: psi_dir(ndmx,2) ! auxiliary function.
  REAL(DP) :: ze2, jnor, thrdum=0.0_dp
  INTEGER  :: n, i, nstop
    =#

    thresh0 = 1.0e-10
    #XXX This is only for rel==1
    # CALL lschps( 3, zed, thrdum, grid, n, 1, l, e, vpot, psi_out, nstop )
    mode = 3
    Nrmesh = grid.Nrmesh
    e, nstop = lschps!( mode, Zval, thresh0, 
        grid, 1, ℓ, E_in, Vpot, psi_out; nin_ = Nrmesh
    )


    # fix arbitrarily the norm at the cut-off radius equal to (about) 0.5**2
    nrm1 = 0.0
    for ir in 1:idx_rc
        nrm1 += grid.dx * grid.r[ir] * psi_out[ir]^2
    end
    nrm1 = sqrt(nrm1)
    for ir in 1:Nrmesh
        psi_out[ir] = psi_out[ir] * 0.5/nrm1
    end
    #IF( rel==2 ) THEN
    #    DO n = 1,grid%mesh
    #        psi_out_rel(n) = psi_out_rel(n)*0.5_dp/jnor
    #    ENDDO
    #ENDIF
    #
    # Set the wavefunction to zero if it diverges too much. In any case
    # this part of the scattering wavefunction is not used and generates 
    # noise in the code.
    for ir in (idx_rc+1):Nrmesh
        if abs(psi_out[ir]) > 1e9
            for irr in ir:Nrmesh
                psi_out[irr] = 0.0
                #IF( rel==2 ) then
                #    psi_out_rel(i) = 0.0_DP
                #endif
            end
        end
    end
    return
end

