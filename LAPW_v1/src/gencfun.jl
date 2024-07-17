function gencfun(pw, atoms, ffacg)
    
    Ng = pw.gvec.Ng
    G = pw.gvec.G
    Natoms = atoms.Natoms
    atpos = atoms.positions
    atm2species = atoms.atm2species

    # allocate global characteristic function arrays
    cfunig = zeros(ComplexF64, Ng)
    
    cfunig[1] = 1.0
    # Loop over all atoms
    for ia in 1:Natoms
        isp = atm2species[ia]
        v1 = atpos[1,ia]
        v2 = atpos[2,ia]
        v3 = atpos[3,ia]
        for ig in 1:Ng
            # structure factor
            t1 = G[1,ig]*v1 + G[2,ig]*v2 + G[3,ig]*v3
            z1 = cos(t1) - im*sin(t1)
            # add to characteristic function in G-space
            cfunig[ig] = cfunig[ig] - ffacg[ig,isp]*z1
        end
    end

    ctmp = zeros(ComplexF64, Npoints)
    for ig in 1:Ng
        ip = pw.gvec.idx_g2r[ig]
        ctmp[ip] = cfunig[ig]
    end
    # Fourier transform to real-space
    G_to_R!(pw, ctmp)
    cfunir .= real.(ctmp)

    return return cfunig, cfunir
end
