
function rhomagk(
    atoms, atsp_vars,
    mt_vars, pw,
    apwalm, evecfv, evecsv,
    rhomt, rhoir
)
# FIXME: need magir and magmt for spinpol
    Natoms = atoms.Natoms
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    npcmt = mt_vars.npcmti
    
    # loop over all atoms
    for ia in 1:Natoms
        isp = atm2species[ia]
        npc = npcmt[isp]
        # de-phasing factor for spin-spirals
        #if ssdph then 
        #    t1 = -0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
        #    zq(1) = cmplx(cos(t1),sin(t1),8)
        #    zq(2) = conjg(zq(1))
        #end
        done[:,:] .= false
        # loop over all second-variational state
        for j in 1:nstsv
            wo = occsvp[j] # occupation number
            # skip this state if it is empty or nearly empty
            if abs(wo) < epsocc
                continue
            end
            #
            wo = wo*wppt
            #
            # if using 2nd variational scheme
            if tevecsv
                # generate spinor wavefunction from second-variational eigenvectors
                i = 0
                for ispn in 1:nspinor
                    jspn = jspnfv[ispn]
                    wfmt3[1:npc,ispn] = 0.0
                    for ist in 1:nstfv
                        i = i + 1
                        z1 = evecsv[i,j]
                        if abs(real(z1)) + abs(imag(z1)) > epsocc
                            #IF(ssdph) z1=z1*zq(ispn)
                            if !done[ist,jspn]
                                wavefmt(lradstp, ias, ngp(jspn), apwalm(:,:,:,ias,jspn), evecfv(:,ist,jspn), wfmt2)
                                # convert to spherical coordinates
                                backward_SHT!(mt_vars, isp, wfmt2, wfmt1)
                                done[ist,jspn] = true
                            end # if
                            # add to spinor wavefunction
                            #CALL zaxpy(npc,z1,wfmt1(:,ist,jspn),1,wfmt3(:,ispn),1)
                            wfmt3[:,ispn] += wfmt1[1:npc,ist,jspn] * wfmt3[1:npc,ispn]
                        end
                    end # do
                end # do
            else 
                # not using 2nd variational scheme
                #
                # spin-unpolarised wavefunction
                wavefmt(lradstp, ias, ngp, apwalm[:,:,:,ias,1], evecfv[:,j,1], wfmt2)
                # The result is stored in wfmt2
                #
                # convert to spherical coordinates
                backward_SHT!(mt_vars, isp, wfmt2, wfmt1)
            end
            #
            # add to density and magnetisation
            if spinpol
                error("Not supported yet")
                # spin-polarised
                if ncmag
                    # non-collinear
                    #rhomagk_rmk1(npc, wo, wfmt3, wfmt3[:,2], rhomt[:,ias], magmt[:,ias,1], magmt[:,ias,2], magmt[:,ias,3])
                else
                    # collinear
                    rhomagk_rmk2(npc, wo, wfmt3, wfmt3[:,2], rhomt[ia], magmt[:,ias,1])
                end
            else
                # spin-unpolarised
                rhomagk_rmk3(npc, wo, wfmt3, rhomt[:,ia])
            end 
  
        end # over states
  
    end # over atoms


    #------------------------------------------------!
    #     interstitial density and magnetisation     !
    #------------------------------------------------!
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    nspinor = elec_chgst.nspinor
    #
    wfir = zeros(ComplexF64, Npoints, nspinor)
    #
    # loop over all states
    #
    for j in 1:nstsv
        wo = occsvp[j]
        if abs(wo) < epsocc
            continue
        end
        wo = wo*wppt/CellVolume
        fill!(wfir, 0.0)
        #
        if tevecsv
            # generate spinor wavefunction from second-variational eigenvectors
            i = 0 # idx for accessing basis function in evecsv
            for ispn in 1:nspinor
                jspn = jspnfv[ispn]
                for ist = 1:nstfv
                    i = i + 1
                    z1 = evecsv[i,j]
                    if abs(real(z1)) + abs(imag(z1)) > epsocc # XXX why check with epsocc?
                        for igp in 1:ngp(jspn)
                            ifg = igfft[igpig[igp,jspn]]
                            wfir[ifg,ispn] += z1*evecfv[igp,ist,jspn]
                        end
                    end # if
                end
            end
        else
            # spin-unpolarised wavefunction
            for igw in 1:Ngwk
                ip = idx_gw2r[igw]
                wfir[ip,1] = evecfv[igp,j,1]
            end
        end
    
        # Fourier transform wavefunction to real-space
        for ispn in 1:nspinor
            #zfftifc(3,ngridg,1,wfir(:,ispn))
            @views G_to_R!(pw, wfir[:,ispn])
        end
        # add to density and magnetisation
        if spinpol 
            # spin-polarised
            if ncmag
                # non-collinear
                @views rhomagk_rmk1!(ngtot, wo, wfir, wfir[:,2], rhoir, magir, magir[:,2], magir[:,3])
            else
                # collinear
                @views rhomagk_rmk2!(ngtot, wo, wfir, wfir[:,2], rhoir, magir)
            end
        else
            # XXX We pass full FFT grid array here, so ngtot -> Npoints
            # spin-unpolarised
            rhomagk_rmk3!(Npoints, wo, wfir, rhoir)
        end
    end # nstsv

    return
end


function rhomagk_rmk1!(n, wo, wf1, wf2, rho, mag1, mag2, mag3)
    for i in 1:n
        z1 = wf1[i]
        z2 = wf2[i]
        t1 = real(z1)^2 + imag(z1)^2
        t2 = real(z2)^2 + imag(z2)^2
        z1 = conj(z1)*z2
        mag1[i] = mag1[i] + wo2*dble(z1)
        mag2[i] = mag2[i] + wo2*aimag(z1)
        mag3[i] = mag3[i] + wo*(t1-t2)
        rho[i] = rho[i] + wo*(t1 + t2)
    end
    return
end


function rhomagk_rmk2!(n, wo, wf1, wf2, rho, mag)
    for i in 1:n
        t1 = real(wf1[i])^2 + imag(wf1[i])^2
        t2 = real(wf2[i])^2 + imag(wf2[i])^2
        mag[i] += wo*(t1 - t2)
        rho[i] += wo*(t1 + t2)
    end
    return
end

# used for non-spin-polarized case
function rhomagk_rmk3!(n, wo, wf, rho)
    for i in 1:n
        rho[i] += wo*( real(wf[i])^2 + imag(wf[i])^2 )
    end
    return
end