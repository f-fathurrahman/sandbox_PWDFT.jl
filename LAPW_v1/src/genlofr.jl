function genlofr!(atoms, eqatoms, mt_vars, apwlo_vars)

    #=
    USE m_atoms, ONLY: natmmax, natoms, nspecies, idxas
  USE m_symmetry, ONLY: eqatoms
  USE m_apwlo, ONLY: nplorb, nlorb, lorbordmax, nlomax, lofr, lorbl, &
               lorbord, deapwlo, lorbdm, lorbe
  USE m_muffin_tins, ONLY: nrmtmax, nrmt, nrmti, rlmt, rmt, lmmaxo, lmmaxi
  USE m_constants, ONLY: y00, solsc
  USE m_density_pot_xc, ONLY: vsmt
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ia,ja,ias,jas
  INTEGER :: nr,nri,ir,i
  INTEGER :: ilo,jlo,io,jo
  INTEGER :: nn,l,info
  REAL(8) :: e,t1
  ! automatic arrays
  LOGICAL :: done(natmmax)
  INTEGER :: ipiv(nplorb)
  REAL(8) :: vr(nrmtmax), fr(nrmtmax)
  REAL(8) :: p0(nrmtmax,lorbordmax), p1(nrmtmax)
  REAL(8) :: q0(nrmtmax), q1(nrmtmax), ep0(nrmtmax,lorbordmax)
  REAL(8) :: p0s(nrmtmax,nlomax), ep0s(nrmtmax,nlomax)
  REAL(8) :: xa(nplorb), ya(nplorb)
  REAL(8) :: a(nplorb,nplorb), b(nplorb)
  ! external functions
  REAL(8) splint, polynm
  external splint, polynm
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
        for ir in (nri+1):nr
            vr[ir] = vsmt[ia][i]*y00
            i = i + lmmaxo
        end
        #
        for ilo in 1:nlorb[isp]
            l = lorbl[isp][ilo]
            for jo in 1:lorbord[isp][ilo]
                # linearisation energy accounting for energy derivative
                E = lorbe[ia][ilo][jo] + lorbdm[isp][ilo][jo]*deapwlo
                # integrate the radial Schrodinger equation
                @views rgrid = rlmt[isp][:,1]
                @views p0view = p0[:,jo]
                nn, E = rschrodint!(l, E, rgrid, vr, p0view, p1, q0, q1)
                ep0[1:nr,jo] .= E*p0[1:nr,jo]
                # normalise radial functions
                fr[1:nr] = p0[1:nr,jo]^2
                t1 = splint(nr, rgrid, fr)
                t1 = 1.0/sqrt(abs(t1))
                #CALL dscal(nr,t1,p0(:,jo),1)
                p0[:,jo] .*= t1
                #CALL dscal(nr,t1,ep0(:,jo),1)
                ep0[:,jo] .*= t1
                # set up the matrix of radial derivatives
                for i in 1:nplorb
                    ir = nr - nplorb + i
                    xa[i] = rlmt[isp][ir,1]
                    ya[i] = p0[ir,jo]*rlmt[isp][ir,-1]
                end
                for io in 1:lorbord[isp][ilo]
                    a[io,jo] = polynm(io-1, nplorb, xa, ya, rmt[isp])
                end 
            end 
            # set up the target vector
            b[:] = 0.0
            b[lorbord[isp][ilo]] = 1.0
            #CALL dgesv(lorbord(ilo,is), 1, a, nplorb, ipiv, b, nplorb, info)
            b[:] = A\b[:]
            #IF(info != 0) goto 10
            # generate linear superposition of radial functions
            p0s[:,ilo] .= 0.0
            ep0s[:,ilo] = 0.0
            for io in 1:lorbord[isp][ilo]
                t1 = b[io]
                #CALL daxpy(nr,t1,p0(:,io),1,p0s(:,ilo),1)
                p0s[:,ilo] .+= t1*p0[:,io]
                #CALL daxpy(nr,t1,ep0(:,io),1,ep0s(:,ilo),1)
                ep0s[:,ilo] .+= t1*ep0[:,io]
            end 
            # normalize radial functions
            fr[1:nr] .= p0s[1:nr,ilo].^2
            t1 = splint(nr, rgrid, fr)
            t1 = 1.0/sqrt(abs(t1))
            #CALL dscal(nr,t1,p0s(:,ilo),1)
            p0s[:,ilo] .*= t1
            #CALL dscal(nr,t1,ep0s(:,ilo),1)
            ep0s[:,ilo] .*= t1
            # subtract linear combination of previous local-orbitals with same l
            for jlo in 1:(ilo-1)
                if lorbl[isp][jlo] == l
                    fr[1:nr] .= p0s[1:nr,ilo] .* p0s[1:nr,jlo]
                    t1 = -splint(nr, rgrid, fr)
                    #CALL daxpy(nr,t1,p0s(:,jlo),1,p0s(:,ilo),1)
                    p0s[:,ilo] .+= t1*p0s[:,jlo]
                    #CALL daxpy(nr,t1,ep0s(:,jlo),1,ep0s(:,ilo),1)
                    ep0s[:,ilo] .+= t1*ep0s[:,ilo]
                end
            end
            # normalize radial functions again
            fr[1:nr] = p0s[1:nr,ilo].^2
            t1 = splint(nr, rgrid, fr)
            t1 = abs(t1)
            if t1 < 1e-25
                error("Degenerate LO radial functions")
            end
            t1 = 1.0/sqrt(t1)
            # CALL dscal(nr,t1,p0s(:,ilo),1)
            p0s[:,ilo] .*= t1
            #CALL dscal(nr,t1,ep0s(:,ilo),1)
            ep0s[:,ilo] .*= t1
            # divide by r and store in global array
            for ir in 1:nr
                t1 = rlmt[isp][ir,-1]
                lofr[ia][ilo][ir,1] = t1*p0s[ir,ilo]
                lofr[ia][ilo][ir,2] = t1*ep0s[ir,ilo]
            end
        end # for ilo
        done[ia] = true
        # copy to equivalent atoms
        for ja in 1:Natoms
            if !done[ja] && eqatoms[ia,ja]
                for ilo in 1:nlorb[isp]
                    #CALL dcopy(nr,lofr(:,1,ilo,ias),1,lofr(:,1,ilo,jas),1)
                    lofr[ja][ilo][:,1] .= lofr[ia][ilo][:,1]
                    #CALL dcopy(nr,lofr(:,2,ilo,ias),1,lofr(:,2,ilo,jas),1)
                    lofr[ja][ilo][:,2] .= lofr[ia][ilo][:,2]
                end 
                done[ja] = true
            end
        end
    end
    return
end
