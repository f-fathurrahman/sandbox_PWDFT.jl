function usnldiag!(
    ik, ispin, atoms, pw, pspotNL,
    H_diag::Vector{Float64}, S_diag::Vector{Float64}
)
    Ngwk = pw.gvecw.Ngw[ik]
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    if any(pspotNL.are_ultrasoft)
        qq_at = pspotNL.qq_at
    else
        nhm = pspotNL.nhm
        qq_at = zeros(Float64,nhm,nhm,Natoms)
    end
    Deeq = pspotNL.Deeq
    betaNL = pspotNL.betaNL
    nh = pspotNL.nh
    indv_ijkb0 = pspotNL.indv_ijkb0

    ps1 = zeros(ComplexF64,2)
    ps2 = zeros(ComplexF64,2)

    for igk in 1:Ngwk
        S_diag[igk] = 1.0
    end

    println("sum Deeq = ", sum(Deeq))
    println("sum qq_at = ", sum(qq_at))
    println("sum betaNL[ik] = ", sum(betaNL[ik]))

    for isp in 1:Nspecies, ia in 1:Natoms

        if atm2species[ia] != isp
            continue
        end
      
        # Loop over projectors
        for ih in 1:nh[isp]
          
            ikb = indv_ijkb0[ia] + ih
          
            ps1[1] = Deeq[ih,ih,ia,ispin]
            ps2[1] = qq_at[ih,ih,ia]

            for igk in 1:Ngwk
                ar = betaNL[ik][igk,ikb]*conj(betaNL[ik][igk,ikb])
                H_diag[igk] = H_diag[igk] + real( ps1[1] * ar )
                S_diag[igk] = S_diag[igk] + real( ps2[1] * ar )
            end
            
            #if upf(nt)%tvanp .or. upf(nt)%is_multiproj ) THEN
            # nh should be larger than 1 ?
            
            if !pspotNL.are_ultrasoft[isp]
                continue
            end

            for jh in 1:nh[isp]
                
                if jh == ih
                    continue
                end
                
                jkb = indv_ijkb0[isp] + jh
                ps1[1] = Deeq[ih,jh,ia,ispin]
                ps2[1] = qq_at[ih,jh,ia]

                for igk in 1:Ngwk
                    ar = betaNL[ik][igk,ikb]*conj(betaNL[ik][igk,ikb])
                    H_diag[igk] = H_diag[igk] + real( ps1[1] * ar )
                    S_diag[igk] = S_diag[igk] + real( ps2[1] * ar )
                end
            end

        end
    
    end

    return
end


# No overlap matrix
function usnldiag!(
    ik, ispin, atoms, pw, pspotNL,
    H_diag::Vector{Float64}
)
    Ngwk = pw.gvecw.Ngw[ik]
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    Deeq = pspotNL.Deeq
    betaNL = pspotNL.betaNL
    nh = pspotNL.nh
    indv_ijkb0 = pspotNL.indv_ijkb0

    ps1 = zeros(ComplexF64,2)

    for isp in 1:Nspecies, ia in 1:Natoms

        if atm2species[ia] != isp
            continue
        end
      
        # Loop over projectors
        for ih in 1:nh[isp]
          
            ikb = indv_ijkb0[ia] + ih
          
            ps1[1] = Deeq[ih,ih,ia,ispin]

            for igk in 1:Ngwk
                ar = betaNL[ik][igk,ikb]*conj(betaNL[ik][igk,ikb])
                H_diag[igk] = H_diag[igk] + real( ps1[1] * ar )
            end
            
            #if upf(nt)%tvanp .or. upf(nt)%is_multiproj ) THEN
            # nh should be larger than 1 ?
            
            if !pspotNL.are_ultrasoft[isp]
                continue
            end

            for jh in 1:nh[isp]
                
                if jh == ih
                    continue
                end
                
                jkb = indv_ijkb0[isp] + jh
                ps1[1] = Deeq[ih,jh,ia,ispin]

                for igk in 1:Ngwk
                    ar = betaNL[ik][igk,ikb]*conj(betaNL[ik][igk,ikb])
                    H_diag[igk] = H_diag[igk] + real( ps1[1] * ar )
                end
            end

        end
    
    end

    return
end
