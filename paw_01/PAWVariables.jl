mutable struct PAWVariables

end

function PAWVariables(paw_data)


end # PAWVariables



# DRAFT
function init_PAW_atomic_becsum!(atoms, pspots, pspotNL; Nspin=1, noise=0.0)

    magn0 = zeros(atoms.Nspecies)

    for ia in 1:Natoms
        #
        isp = atm2species[ia]
        psp = pspots[isp]
        #
        if pspots[isp].is_paw
            ijh = 1
            for ih in 1:nh[isp]
                #
                nb = indv[ih,isp]
                occnb = psp.paw.oc[nb]
                #
                if Nspin == 1
                    becsum[ijh,ia,1] = occnb / ( 2*nhtol[ih,isp] + 1 )
                else # Nspin == 2
                    becsum[ijh,ia,1] = 0.5*(1.0 + magn0[isp] ) * occnb / ( 2*nhtol[ih,isp] + 1 )
                    becsum[ijh,ia,2] = 0.5*(1.0 - magn0[isp] ) * occnb / ( 2*nhtol[ih,isp] + 1 )
                end
                ijh = ijh + 1
                for jh in (ih + 1):nh[isp]
                    for ispin in 1:Nspin
                        if noise > 0.0
                            becsum[ijh,ia,ispin] += noise *2.0*( 0.5 - rand() )
                        end
                    end
                    ijh = ijh + 1
                end # jh
            end # ih
        end # if paw
    end # ia
    
    # copy becsum in scf structure and symmetrize it
    #rho%bec(:,:,:) = becsum(:,:,:)
    

    #CALL PAW_symmetrize( rho%bec )

    return
end