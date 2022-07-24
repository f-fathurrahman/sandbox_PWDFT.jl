function _calc_becsum( ik, ispin, Ham, psiks )

    Natoms = Ham.atoms.Natoms
    Nspecies = Ham.atoms.Nspecies
    atm2species = Ham.atoms.atm2species

    Nspin = Ham.electrons.Nspin
    Focc = Ham.electrons.Focc
    
    nhm = Ham.pspotNL.nhm
    nh = Ham.pspotNL.nh
    indv_ijkb0 = Ham.pspotNL.indv_ijkb0

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    ikspin = ik + (ispin-1)*Nkpt
    psi = psiks[ikspin]

    Nstates = size(psi, 2)

    betaNL_psi = Ham.pspotNL.betaNL[ik]' * psi

    Nbecsum = Int64( nhm * (nhm + 1)/2 )
    becsum = zeros(Float64, Nbecsum, Natoms, Nspin)

    for isp in 1:Nspecies

        # Skip if this species does not use ultrasoft
        if !Ham.pspots[isp].is_ultrasoft
            continue
        end
        
        # These are used for GEMM.
        # They can be removed
        auxk1 = zeros(ComplexF64, Nstates, nh[isp])
        auxk2 = zeros(ComplexF64, Nstates, nh[isp])
        #
        # In becp=<vkb_i|psi_j> terms corresponding to atom ia of type isp
        # run from index i=indv_ijkb0[ia]+1 to i=indv_ijkb0[ia] + nh[isp]
        for ia in 1:Natoms
            
            if atm2species[ia] != isp
                continue
            end

            # sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
            # copy into aux1, aux2 the needed data to perform a GEMM
            for ih in 1:nh[isp]
                ikb = indv_ijkb0[ia] + ih
                for ist in 1:Nstates
                    auxk1[ist,ih] = betaNL_psi[ikb,ist]
                    auxk2[ist,ih] = Focc[ist,ikspin]*betaNL_psi[ikb,ist]
                end
            end
            #
            # only the real part is computed
            #
            #CALL DGEMM ( 'C', 'N', nh(np), nh(np), 2*this_bgrp_nbnd, &
            #     1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
            #     0.0_dp, aux_gk, nh(np) )
            aux_gk = real(auxk1' * auxk2)

            # TODO: calculated aux_gk directly without auxk1 and auxk2
            
            # copy output from GEMM into desired format
            #
            ijh = 0
            for ih in 1:nh[isp], jh in ih:nh[isp]
                ijh = ijh + 1
                # nondiagonal terms summed and collapsed into a
                # single index (matrix is symmetric wrt (ih,jh))
                if jh == ih
                    becsum[ijh,ia,ispin] += aux_gk[ih,jh]
                else
                    becsum[ijh,ia,ispin] += aux_gk[ih,jh]*2.0
                end
            end
        end
    end

    return becsum

end
