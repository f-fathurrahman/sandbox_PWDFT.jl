function olprad!( atoms, mt_vars, apwlo_vars, oalo, ololo )
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

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    nrmt = mt_vars.nrmt
    nrmtmax = maximum(nrmt)
    wrmt = mt_vars.wrmt

    nlorb = apwlo_vars.nlorb
    lorbl = apwlo_vars.lorbl
    apword = apwlo_vars.apword
    apwfr = apwlo_vars.apwfr
    lofr = apwlo_vars.lofr

    fr = zeros(Float64, nrmtmax)

    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        #-------------------------------------
        #     APW-local-orbital integrals     
        #-------------------------------------
        for ilo in 1:nlorb[isp]
            l = lorbl[isp][ilo]
            for io in 1:apword[isp][l]
                fr[1:nr] .= apwfr[ia][l][io][1:nr,1] .* lofr[ia][ilo][1:nr,1]
                oalo[ia][io,ilo] = dot( wrmt[isp][1:nr], fr[1:nr] )
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
                    ololo[ia][ilo,jlo] = dot( wrmt[isp][1:nr], fr[1:nr] )
                end 
            end
        end
    end
    return
end
