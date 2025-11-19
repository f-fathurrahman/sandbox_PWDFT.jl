function rhomagsh!(atoms, mt_vars, rhomt; magmt=nothing)

    if isnothing(magmt)
        ndmag = 0
    else
        ndmag = size(magmt[1], 2)
    end

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    npcmt = mt_vars.npcmt
    npcmtmax = maximum(npcmt)

    rfmt = zeros(Float64, npcmtmax)
    for ia in 1:Natoms
        isp = atm2species[ia]
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        npc = npcmt[isp]
        # convert the density to spherical harmonics
        #println("size rhomt[ia] = ", size(rhomt[ia]))
        #println("npmct = ", npcmt)
        rfmt[1:npc] = rhomt[ia][1:npc]
        forward_SHT!( mt_vars, isp, rfmt, rhomt[ia], coarse=true )
        # convert magnetisation to spherical harmonics
        for idm in 1:ndmag
            rfmt[1:npc] = magmt[ia][1:npc]
            @views forward_SHT!( mt_vars, isp, rfmt, magmt[ia][:,idm], coarse=true )
        end
    end
    return
end

