
function calc_forces_nlcc!( F_nlcc )
    # Calculates the NLCC contribution to the force

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies

    Npoints = prod(pw.Ns)

    # early return if there is no pseudopot with nlcc  
    IF( ANY( upf(1:ntyp)%nlcc ) ) GOTO 15

    fact = pw.CellVolume # 2*pw.CellVolume if using gamma only
  
    #
    # recalculate the exchange-correlation potential
    #
    Vxc = zeros(Float64, Npoints, Nspin)
    # CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )

    psic = zeros(ComplexF64, Npoints)
    @assert Nspin == 1
    for ir in 1:Npoints
        psic[ir] = Vxc[ir,1]
    end

    #fwfft('Rho', psic, dfftp )
    #
    # psic contains now Vxc(G)
    #
    rhocg = zeros(Float64, Ngl)
    #
    # core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
    # G = 0 term gives no contribution
    #
    for isp in 1:Nspecies
        
        # skip if no nlcc
        if !psp.is_nlcc
            continue
        end

        #CALL drhoc( ngl, gl, omega, tpiba2, msh(nt), rgrid(nt)%r, &
        #          rgrid(nt)%rab, upf(nt)%rho_atc, rhocg )
        
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            for ig in 2:Ng
                igl = idx_g2shells[ig]
                ip = idx_g2r[ig]
                GX = G[1,ig]*X[1,ia] + G[2,ig]*X[2,ia] + G[3,ig]*X[3,ia]
                Sf = sin(GX) + im*cos(GX)
                F_nlcc[:,ia] .+= fact * G[:,ig] * rhocg[igl] * conj(psic[ip]) * Sf
            end
        end
    end
    return
end
