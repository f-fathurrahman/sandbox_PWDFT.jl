function olprad!()
#=
  USE m_atoms, ONLY: natmtot, idxis
  USE m_muffin_tins, ONLY: nrmtmax, nrmt, wrmt
  USE m_apwlo, ONLY: lofr, nlorb, lorbl, apword, apwfr
  USE m_hamiltonian, ONLY: ololo, oalo
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ias,nr
  INTEGER :: ilo,jlo,l,io
  ! automatic arrays
  REAL(8) :: fr(nrmtmax)
=#
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        #-------------------------------------
        #     APW-local-orbital integrals     
        #-------------------------------------
        for ilo in 1:nlorb[isp]
            l = lorbl[isp][ilo]
            for io in 1:apword[isp][l]
                fr[1:nr] .= apwfr[ia][l][ip][1:nr,1] .* lofr[ia][ilo][1:nr,1]
                oalo[ia][io,ilo] = dot( wrmt[1:nr,isp], fr[1:nr] )
            end
        end
        #----------------------------------------------
        #     local-orbital-local-orbital integrals
        #-----------------------------------------------
        for ilo in 1:nlorb[isp]
            l = lorbl[isp][ilo]
            for jlo in 1:nlorb[isp]
                if lorbl[isp][jlo] == l
                    fr[1:nr] .= lofr[ia][ilo][1:nr,1] .* lofr[ia][jlo][1:nr,1]
                    ololo[ia][ilo,jlo] = dot( wrmt[1:nr,is], fr[1:nr] )
                end 
            end
        end
    end
    return
end
