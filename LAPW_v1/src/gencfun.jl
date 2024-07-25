function gencfun(pw, atoms, ffacg; gvec_full=nothing)

    if isnothing(gvec_full)
        gvec = pw.gvec
    else
        gvec = gvec_full
    end
    
    Ng = gvec.Ng
    G = gvec.G

    # This should the same as prod(gvec_full.Ng) in case gvec_full is used
    Npoints = prod( Ns)

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
            # XXX: Elk uses more G-vectors for cfunig and ffacg
            # structure factor
            t1 = G[1,ig]*v1 + G[2,ig]*v2 + G[3,ig]*v3
            z1 = cos(t1) - im*sin(t1)
            # add to characteristic function in G-space
            cfunig[ig] = cfunig[ig] - ffacg[ig,isp]*z1
        end
    end

    ctmp = zeros(ComplexF64, Npoints)
    for ig in 1:Ng
        ip = gvec.idx_g2r[ig]
        ctmp[ip] = cfunig[ig]
    end
    # Fourier transform to real-space
    G_to_R!(pw, ctmp)
    cfunir = real.(ctmp)*Npoints # scale

    return return cfunig, cfunir
end
