function genbs!(
    atoms, mt_vars,
    cfunir,
    bfcmt, bfieldc,
    bxcmt, bxcir,
    bsmt, bsir
)
    #=
  ! local variables
  INTEGER :: idm,is,ia,ias
  INTEGER :: nrc,nrci,npc
  REAL(8) :: cb,t1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:)
    =#

    # check spinpol
    if ndmag == 0
        return
    end
    
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    npcmt = mt_vars.npcmt

    solsc = 137.035999084 # XXX: not scaled
    gfacte = 2.00231930436256
    # coupling constant of the external field (g_e/4c)
    cb = gfacte/(4.0*solsc)

    # muffin-tin Kohn-Sham field     !
    rfmt = zeros(Float64, npcmtmax)
    
    for ia in 1:Natoms
        isp = atm2species[ia]
        nrc = nrcmt[ip]
        nrci = nrcmti[isp]
        npc = npcmt[isp]
        # exchange-correlation magnetic field in spherical coordinates
        for idm in 1:ndmag
            @views rf_mt_f_to_c( mt_vars, nrc, nrci, bxcmt[ia][:,idm], rfmt )
            # rbsht(nrc,nrci,rfmt,bsmt(:,ias,idm))
            @views backward_SHT!( mt_vars, isp, rfmt[ia], bsmt[ia][:,idm] )
        end 
        # add the external magnetic field
        t1 = cb*( bfcmt[3,ia] + bfieldc[3] )
        bsmt[ia][1:npc,ndmag] = bsmt[ia][1:npc,ndmag] + t1
        if ncmag 
            for idm in 1:2
                t1 = cb*( bfcmt[idm,ia] + bfieldc[idm] )
                bsmt[ia][1:npc,idm] = bsmt[ia][1:npc,idm] + t1
            end
        end # if
    end # do
    # rfmt is not used anymore

    # interstitial Kohn-Sham magnetic field
    for idm in 1:ndmag
        if ncmag
            t1 = cb*bfieldc[idm]
        else
            t1 = cb*bfieldc[3]
        end
        # multiply by characteristic function
        bsir[:,idm] = ( bxcir[:,idm] + t1)*cfunir[:]
    end
    
    # add the magnetic dipole field if required
    #if tbdip
    #    # bdipole()
    #end

    return
end
