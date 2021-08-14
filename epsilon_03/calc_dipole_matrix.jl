function calc_dipole_matrix!(
    pw, psi,
    Nspin, Focc, ik, M; metal_like=false
)

    @assert Nspin == 1

    fill!( M, 0.0 + im*0.0 )

    Nstates = size(psi,2)
    Npw = size(psi,1)

    # XXX Access at ik=1 (memory-saving version)
    idx_gw2g = pw.gvecw.idx_gw2g[1]
    G = pw.gvec.G
    @views k = pw.gvecw.kpoints.k[:,1]

    FULL_OCC = 2.0 # for non-spinpol case

    # XXX Focc is only needed to determine which states we should calculate

    for jst in 1:Nstates
        if Focc[jst,1] >= FULL_OCC # skip if states are occupied
            continue
        end
        for ist in 1:Nstates
            #
            if ist == jst
                continue
            end
            #
            if Focc[ist,1] >= 0.5e-4*FULL_OCC # occupied states
                for igw in 1:Npw
                    ig = idx_gw2g[igw]
                    caux = conj(psi[igw,ist])*psi[igw,jst]
                    for i in 1:3
                        M[i,ist,jst] = M[i,ist,jst] + ( G[i,ig] + k[i] ) * caux
                    end
                end
            end
        end
    end
  
    # The diagonal terms are taken into account only if the system is treated like a metal,
    # not in the intraband term. Because of this we can recalculate the diagonal
    # component of the dipole
    # tensor directly as we need it for the intraband term, without interference with interband one.
    if metal_like
        for ist in 1:Nstates
            for igw in 1:Npw
                ig = idx_gw2g[igw]
                caux = conj(psi[igw,ist])*psi[igw,ist]
                for i in 1:3
                    M[i,ist,ist] = M[i,ist,ist] + ( G[i,ig] + k[i] ) * caux
                end
            end
        end
    end

    return
end