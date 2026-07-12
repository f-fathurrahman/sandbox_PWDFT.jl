function op_Vexx!(Ham, exx, psi, Vpsi)

    Npoints_exx = prod(exx.Ns)

    temppsic = zeros(ComplexF64, Npoints_exx)
    res = zeros(ComplexF64, Npoints_exx)

    fac = zeros(Float64, Npoints_exx)
    rhoc = zeros(ComplexF64, Npoints_exx)
    vc = zeros(ComplexF64, Npoints_exx)
    
    current_ik = Ham.ik # this is current index of the input wavefunction
    xkp = Ham.pw.gvecw.kpoints.k[:,current_ik]
    exx_alpha = 1.0 # FIXME
    
    Ncols = size(psi, 2)
    for ic in 1:Ncols
        fill!(temppsic, 0.0)
        # Bring psi(:,im) to real space using invfft, using exx grid
        for igw in 1:Ngwk
            ip = idx_gw2r[ik][igw]
            temppsic[ip] = psi[ig,ic]
        end
        do_fft!(planbw, Ns, temppsic)
        #
        fill!(res, 0.0)
        #
        for iq in 1:nqs
            ikq = index_xkq[current_ik,iq]
            ik = index_xk[ikq]
            xkq = xkq_collect[:,ikq]
            # calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
            g2_convolution!( pw, exx, gt, xkp, xkq, iq, current_ik )
            fill!(fac, 0.0)
            for ig in 1:Ng
                ip = idx_g2r[ig]
                fac[ip] = coulomb_fac[ig, iq, current_ik]
            end
            for ist in 1:Nstates
                if abs(x_occupation[ist,ik]) < eps_occ
                    continue
                end
                for ip in 1:Npoints_exx
                    rhoc[ip] = conj(exx.buff[ip,ibnd,ikq])*temppsic[ip] / CellVolume
                end
                # to G-space
                do_fft!(planfw, Ns, rhoc)
                # charge done
                fill!(vc, 0.0)
                #
                # multiply point by points?
                for ip in 1:Npoints_exx
                    vc[ip] = fac[ip] * rhoc[ip]*x_occupation[ibnd,ik]/nqs
                end
                #
                # brings back v in real space
                do_fft!(planbw, Ns, vc)
                #
                # accumulates over bands and k points
                for ip in 1:Npoints_exx
                    result[ip] += vc[ip]*exx.buff[ip,ibnd,ikq]
                end
            end
        end
        #
        # brings back result in G-space
        do_fft!(planfw, Ns, result)
        # adds it to hpsi
        for igw in 1:Ngwk
            ip = idx_gw2r[current_ik][igw]
            Vpsi[igw,ic] -= exx_alpha * result[ip]
        end
    end
    return
end
