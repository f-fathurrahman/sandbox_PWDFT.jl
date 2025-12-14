function init_tab_at(psp, pw)
    #=
    This routine computes a table with the radial Fourier transform 
    of the atomic wavefunctions.

  USE kinds,        ONLY : DP
  USE atom,         ONLY : rgrid, msh
  USE constants,    ONLY : fpi
  USE cell_base,    ONLY : omega
  USE ions_base,    ONLY : ntyp => nsp
  USE us,           ONLY : tab_at, nqx, dq
  USE uspp_param,   ONLY : upf
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum

  IMPLICIT NONE
  INTEGER :: nt, nb, iq, ir, l, startq, lastq, ndm

    REAL(DP), ALLOCATABLE :: aux(:), vchi(:)
  REAL(DP) :: vqint, pref, q
    =#

    ecutwfc = pw.ecutwfc
    CellVolume = pw.CellVolume
    cell_factor = 1.0 # XXX HARDCODED
    dq = 0.01 # XXX HARDCODED
    ndm = psp.Nr_rcut
    aux = zeros(Float64, ndm)
    vchi = zeros(Float64, ndm)
    
    # chiq = radial fourier transform of atomic orbitals chi
    pref = 4π / sqrt(CellVolume)
    # needed to normalize atomic wfcs (not a bad idea in general and 
    # necessary to compute correctly lda+U projections)
    
    Nq = floor( Int64, (sqrt(2*ecutwfc)/dq + 4)*cell_factor )
    Nwfc = psp.Nchi
    tab_at = zeros(Float64, Nq, Nwfc)
    for iwf in 1:Nwfc
        if psp.occ_chi[iwf] >= 0.0
            l = psp.lchi[iwf]
            for iq in 1:Nq
                q = dq * (iq - 1)
                PWDFT.qe_sph_bes!(l, q, psp.r[1:ndm], aux)
                for ir in 1:ndm
                    vchi[ir] = psp.chi[ir,iwf] * aux[ir] * psp.r[ir]
                end
                vqint = PWDFT.integ_simpson( ndm, vchi, psp.rab )
                tab_at[iq,iwf] = vqint * pref
            end
        end # if
    end # iwf
    return tab_at
end
