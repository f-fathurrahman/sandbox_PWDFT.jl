function op_Vexx!(Ham, exx, psi, Vpsi)

    Npoints_exx = prod(exx.Ns)

    temppsic = zeros(ComplexF64, Npoints_exx)
    res = zeros(ComplexF64, Npoints_exx)

    Ng_exx = exx.gvec.Ng
    fac = zeros(Float64, Ng_exx)
    rhoc = zeros(ComplexF64, Npoints_exx)
    vc = zeros(ComplexF64, Npoints_exx)
    
    current_ik = Ham.ik # this is current index of the input wavefunction
    xkp = Ham.pw.gvecw.kpoints.k[:,current_ik]
    exx_alpha = exx.exx_alpha

    pw = Ham.pw
    xk = pw.gvecw.kpoints.k[:,current_ik]
    Nstates = Ham.electrons.Nstates
    x_occupation = exx.x_occupation
    
    Ncols = size(psi, 2)
    for ic in 1:Ncols
        fill!(temppsic, 0.0)
        # Bring psi(:,im) to real space using invfft, using exx grid
        for igw in 1:pw.gvecw.Ngw[current_ik]
            ip = exx.idx_gw2r[current_ik][igw]
            temppsic[ip] = psi[igw,ic]
        end
        do_fft!(exx.planbw, exx.Ns, temppsic)
        #
        fill!(res, 0.0)
        #
        for iq in 1:exx.nqs
            ikq = exx.index_xkq[current_ik,iq]
            ik = exx.index_xk[ikq]
            xkq = exx.xkq[:,ikq]
            # calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
            fill!(fac, 0.0) # need this?
            g2_convolution!( exx, pw.LatVecs, xk, xkq, fac )
            for ist in 1:Nstates
                if abs(x_occupation[ist,ik]) < exx.eps_occ
                    continue
                end
                for ip in 1:Npoints_exx
                    rhoc[ip] = conj(exx.buff[ip,ist,ikq])*temppsic[ip] / pw.CellVolume
                end
                # to G-space
                do_fft!(exx.planfw, exx.Ns, rhoc)
                # charge done
                fill!(vc, 0.0)
                #
                # multiply point by points?
                for ig in 1:Ng_exx
                    ip = exx.gvec.idx_g2r[ig]
                    vc[ip] = fac[ig] * rhoc[ip]*x_occupation[ist,ik]/exx.nqs
                end
                #
                # brings back v in real space
                do_fft!(exx.planbw, exx.Ns, vc)
                #
                # accumulates over bands and k points
                for ip in 1:Npoints_exx
                    res[ip] += vc[ip]*exx.buff[ip,ist,ikq]
                end
            end
        end
        #
        # brings back result in G-space
        do_fft!(exx.planfw, exx.Ns, res)
        # adds it to hpsi
        for igw in 1:pw.gvecw.Ngw[current_ik]
            ip = exx.idx_gw2r[current_ik][igw] # we need to access exx grid
            Vpsi[igw,ic] -= exx_alpha * res[ip] # the sign is negative
        end
    end
    return
end
