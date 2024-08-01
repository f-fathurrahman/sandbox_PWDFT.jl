function genapwfr!(atoms, mt_vars, apwlo_vars)
    #=
  USE m_atoms, ONLY: natmmax, natoms, idxas, nspecies
  USE m_symmetry, ONLY: eqatoms
  USE m_muffin_tins, ONLY: rlmt, nrmt, nrmti, nrmtmax, lmmaxi, lmmaxo, lmaxapw
  USE m_apwlo, ONLY: apword, apwfr, apwordmax, apwdfr, deapwlo, apwdm, apwe
  USE m_density_pot_xc, ONLY: vsmt
  USE m_constants, ONLY: solsc, y00
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ia,ja,ias,jas
  INTEGER :: nr,nri,ir,i
  INTEGER :: nn,l,io,jo
  REAL(8) :: e,t1
  ! automatic arrays
  LOGICAL :: done(natmmax)
  REAL(8) :: vr(nrmtmax),fr(nrmtmax)
  REAL(8) :: p0(nrmtmax,apwordmax),p1(nrmtmax),p1s(apwordmax)
  REAL(8) :: q0(nrmtmax),q1(nrmtmax),ep0(nrmtmax,apwordmax)
  ! external functions
  REAL(8) splint
  EXTERNAL splint
    =#

    y00 = 0.5/sqrt(pi)

    done = zeros(Bool, Natoms)

    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        if done[ia]
            continue
        end 
        # use spherical part of potential
        i = 1
        for ir in 1:nri
            vr[ir] = vsmt[ia][i]*y00
            i = i + lmmaxi
        end
        #
        for ir in (nri+1):nr
            vr[ir] = vsmt[ia][i]*y00
            i = i + lmmaxo
        end
        #
        for l in 0:lmaxapw
            for io in 1:apword[isp][l]
                # linearisation energy accounting for energy derivative
                # XXX FIXME apwdm indexing apwdm[is][l][io]
                E = apwe[ia][l][io] + apwdm[is][io,l]*deapwlo
                # integrate the radial Schrodinger equation
                @views rgrid = rlmt[isp][:,1]
                @views p0view = p0[:,io]
                nn, E = rschrodint!(l, E, nr, grid, vr, p0view, p1, q0, q1)
                # multiply by the linearisation energy
                ep0[1:nr,io] .= E*p0[1:nr,io]
                # normalize radial functions
                fr[1:nr] = p0[1:nr,io]^2
                t1 = splint(nr, rgrid, fr)
                t1 = 1.0/sqrt(abs(t1))
                #CALL dscal(nr,t1,p0(:,io),1)
                p0[:,io] .*= t1
                #
                p1s[io] = t1*p1[nr]
                #CALL dscal(nr,t1,ep0(:,io),1)
                ep0[:,io] .*= t1
                # subtract linear combination of previous vectors
                for jo in 1:(io-1)
                    fr[1:nr] .= p0[1:nr,io] .* p0[1:nr,jo]
                    t1 = -splint(nr, rgrid, fr)
                    #
                    #CALL daxpy(nr,t1,p0(:,jo),1,p0(:,io),1)
                    p0[:,io] .+= t1*p0[:,jo]
                    #
                    p1s[io] = p1s[io] + t1*p1s[jo]
                    #
                    #CALL daxpy(nr, t1, ep0(:,jo), 1, ep0(:,io), 1)
                    ep0[:,io] .+= t1*ep0[:,jo]
                end 
                # normalize radial functions again
                fr[1:nr] .= p0[1:nr,io].^2
                t1 = splint(nr, rgrid, fr)
                t1 = abs(t1)
                if t1 < 1e-25
                    error("Degenerate APW radial functions") 
                end
                t1 = 1.0/sqrt(t1)
                #CALL dscal(nr,t1,p0(:,io),1)
                p0[:,io] .*= t1
                #
                p1s[io] = t1*p1s[io]
                #
                #CALL dscal(nr,t1,ep0(:,io),1)
                ep0[:,io] .*= t1
                #
                # divide by r and store in global array
                for ir in 1:nr
                    t1 = rlmt[isp][ir,-1]
                    apwfr[ia][l][io][ir,1] = t1*p0[ir,io]
                    apwfr[ia][l][io][ir,2] = t1*ep0[ir,io]
                end
                # derivative at the muffin-tin surface
                apwdfr[ia][l][io] = ( p1s[io] - p0[nr,io]*t1 )*t1
            end
        end
        done[ia] = true
        # copy to equivalent atoms
        for ja in 1:Natoms
            if !done[ja] && eqatoms[ia,ja]
                for l in 0:lmaxapw
                    for io in 1:apword[is][l]
                        apwfr[ja][l][io][:,1] .= apwfr[ia][l][io][:,1]
                        apwfr[ja][l][io][:,2] .= apwfr[ia][l][io][:,2]
                        apwdfr[ja][l][io] = apwdfr[ia][l][io]
                    end 
                end
                done[ja] = true
            end
        end
    end

    return

end
