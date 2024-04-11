function zpotcoul!(
    atoms, atsp_vars,
    mt_vars, pw, 
    zrhoir, # input
    zvclmt, # input output
    zvclir 
)

    #println("enter zpotcoul")

    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    rmt = mt_vars.rmt
    gmaxvr = sqrt(2*pw.ecutrho)
    lmaxo = mt_vars.lmaxo
    Npoints = prod(pw.Ns)

    # Poisson solver pseudocharge density constant
    if Nspecies > 0
        t1 = 0.250*gmaxvr*maximum(rmt)
    else
        t1 = 0.25*gmaxvr*2.0  # FIXME: for which case is this useful?
    end
    #println("t1 = ", t1)
    npsd = max(round(Int64,t1),1)
    lnpsd = lmaxo + npsd + 1

    #println("npsd = ", npsd)
    #println("lnpsd = ", lnpsd)

    # compute (R_mt)^l
    rmtl = OffsetArray( zeros(Float64, lmaxo+4, Nspecies),
        0:lmaxo+3, 1:Nspecies )
    for isp in 1:Nspecies
        rmtl[0,isp] = 1.0
        for l in 1:lmaxo+3
            rmtl[l,isp] = rmtl[l-1,isp]*rmt[isp]
        end
    end

    #println("rmtl = ")
    #display(rmtl); println()

    # compute the multipole moments from the muffin-tin potentials
    npmt = mt_vars.npmt
    npmti = mt_vars.npmti
    lmmaxo = mt_vars.lmmaxo

    qlm = zeros(ComplexF64, lmmaxo, Natoms)
    t1 = 1.0/(4*pi)
    for ia in 1:Natoms
        isp = atm2species[ia]
        i = npmt[isp] - lmmaxo # FIXME: idx for accessing zvclmt
        lm = 0
        for l in 0:lmaxo
            t2 = t1*(2*l+1)*rmtl[l+1,isp]
            for m in -l:l
                lm = lm + 1
                i = i + 1
                qlm[lm,ia] = t2*zvclmt[ia][i]
            end
        end
    end

    #println("qlm = ")
    #display(qlm); println()


    # Fourier transform density to G-space and store in zvclir
    zvclir[:] .= zrhoir[:]
    #println("Before R_to_G! sum zvclir = ", sum(zvclir))
    R_to_G!(pw, zvclir)
    # FIXME: scale zvclir with 1/Npoints?
    zvclir[:] = zvclir[:]/Npoints
    #println("After R_to_G! sum zvclir = ", sum(zvclir))


    # FIXME: need to precompute this?
    G2 = pw.gvec.G2
    y00 = 0.28209479177387814347
    epslat = 1e-6
    Ng = pw.gvec.Ng


    ylmg = zeros(ComplexF64, lmmaxo, Ng)
    genylmg!(mt_vars.lmaxo, pw.gvec.G, ylmg)
    #
    sfacg = calc_sfacg(atoms, pw)

    jlgprmt = OffsetArray(
        zeros(ComplexF64, lnpsd+1, Ng, Nspecies),
        0:lnpsd, 1:Ng, 1:Nspecies )
    genjlgprmt!(rmt, lnpsd, G2, jlgprmt)

    # subtract the multipole moments of the interstitial charge density
    for ia in 1:Natoms
        isp = atm2species[ia]
        #for l in 0,lmaxo
        #    zl[l] = 4*pi*zil[l]*rmtl[l+2,isp]
        #end
        for ig in 1:Ng
            ip = pw.gvec.idx_g2r[ig]
            if sqrt(G2[ig]) > epslat
                z1 = zvclir[ip]*sfacg[ig,ia]/sqrt(G2[ig])
                lm = 0
                for l in 0:lmaxo
                    z2 = jlgprmt[l+1,ig,isp]*z1 * 4*pi* (im)^l *rmtl[l+2,isp]
                    for m in -l:l
                        lm = lm + 1
                        qlm[lm,ia] = qlm[lm,ia] - z2*conj(ylmg[lm,ig])
                    end
                end
            else
                t1 = (4*pi/3)*rmtl[3,isp]*y00
                qlm[1,ia] = qlm[1,ia] - t1*zvclir[ip]
            end
        end
    end

    #println("qlm after subtracting multipole moments of zvclir = ")
    #display(qlm); println()


    # find the smooth pseudocharge within the muffin-tin whose multipoles are the
    # difference between the real muffin-tin and interstitial multipoles
    zlm = zeros(ComplexF64, lmmaxo)
    CellVolume = pw.CellVolume
    t1 = (4*pi/CellVolume)*factnm(2*lnpsd+1,2)
    for ia in 1:Natoms
        isp = atm2species[ia]
        lm = 0
        for l in 0:lmaxo
            t2 = t1/( factnm(2*l+1,2)*rmtl[l,isp] )
            z1 = t2*(-1.0*im)^l # TODO: check zilc
            for m in -l:l
                lm = lm + 1
                zlm[lm] = z1*qlm[lm,ia]
            end
        end
        # add the pseudocharge and real interstitial densities in G-space
        for ig in 1:Ng
            ip = pw.gvec.idx_g2r[ig]
            Glen = sqrt(G2[ig])
            if Glen > epslat
                t2 = Glen*rmt[isp]
                t3 = 1.0/t2^lnpsd
                z1 = t3*zlm[1]*ylmg[1,ig]
                lm = 1
                for l in 1:lmaxo
                    lm = lm + 1
                    z2 = zlm[lm]*ylmg[lm,ig]
                    for m in 1-l:l
                        lm = lm + 1
                        z2 = z2 + zlm[lm]*ylmg[lm,ig]
                    end
                    t3 = t3*t2
                    z1 = z1 + t3*z2
                end
                z2 = jlgprmt[lnpsd,ig,isp]*conj(sfacg[ig,ia])
                zvclir[ip] = zvclir[ip] + z1*z2
            else
                t2 = y00/factnm(2*lnpsd+1,2)
                zvclir[ip] = zvclir[ip] + t2*zlm[1]
            end
        end # loop over Ng
    end
    #println("sum zvclir after smooth multipole = ", sum(zvclir))



    # solve Poisson's equation in G+p-space for the pseudocharge
    zvclir[1] = 0.0 + im*0.0 # ip=1 corresponds to ig=1 or G=0
    for ig in 2:Ng
        ip = pw.gvec.idx_g2r[ig]
        zvclir[ip] = zvclir[ip]*(4π/G2[ig])
    end
    #println("sum zvclir after Poisson solve: ", sum(zvclir))


    #
    # match potentials at muffin-tin boundary by adding homogeneous solution
    #
    npmtmax = maximum(npmt)
    zhmt = zeros(ComplexF64,npmtmax)
    lmaxi = mt_vars.lmaxi
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    rlmt = mt_vars.rlmt
    lmmaxi = mt_vars.lmmaxi
    #
    for ia in 1:Natoms
        isp = atm2species[ia]
        # find the spherical harmonic expansion of the interstitial potential at the
        # muffin-tin radius
        fill!(zlm, 0.0)
        for ig in 1:Ng
            ip = pw.gvec.idx_g2r[ig]
            z1 = 4*pi*zvclir[ip]*sfacg[ig,ia]
            lm = 0
            for l in 0:lmaxo
                z2 = jlgprmt[l,ig,isp]*z1*(im)^l
                for m in -l:l
                    lm = lm + 1
                    zlm[lm] = zlm[lm] + z2*conj(ylmg[lm,ig])
                end
            end
        end
        # calculate the homogenous solution
        i = npmt[isp] - lmmaxo
        lm = 0
        #
        for l in 0:lmaxi
            t1 = 1.0/rmtl[l,isp]
            for m in -l:l
                lm = lm + 1
                i = i + 1
                z1 = t1*( zlm[lm] - zvclmt[ia][i] )
                j = lm
                for ir in 1:nrmti[isp]
                    zhmt[j] = z1*rlmt[isp][ir,l]
                    j = j + lmmaxi
                end 
                for ir in (nrmti[isp]+1):nrmt[isp]
                    zhmt[j] = z1*rlmt[isp][ir,l]
                    j = j + lmmaxo
                end
            end
        end
        #
        for l in (lmaxi+1):lmaxo
            t1 = 1.0/rmtl[l,isp]
            for m in -l:l
                lm = lm + 1
                i = i + 1
                z1 = t1*( zlm[lm] - zvclmt[ia][i] )
                j = npmti[isp] + lm
                for ir in (nrmti[isp]+1):nrmt[isp]
                    zhmt[j] = z1*rlmt[isp][ir,l]
                    j = j + lmmaxo
                end
            end
        end

        @views zvclmt[ia][1:npmt[isp]] = zvclmt[ia][1:npmt[isp]] + zhmt[1:npmt[isp]]
        ss = sum(zvclmt[ia])
        #@printf("sum zvclmt[ia] = (%18.10e,%18.10e)\n", real(ss), imag(ss))
    end


    # Convert zvclir to real space
    G_to_R!(pw, zvclir)
    zvclir[:] = zvclir[:]*Npoints  # FIXME: remove Npoints factor (?)

    ss = sum(abs.(zvclir))
    #@printf("sum abs zvclir after G_to_R: (%18.10f,%18.10f)\n", real(ss), imag(ss))

    return
end
