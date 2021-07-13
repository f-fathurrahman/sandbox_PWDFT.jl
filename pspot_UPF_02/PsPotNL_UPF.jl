#mutable struct PsPotNL_UPF
#end

function PsPotNL_UPF( atoms::Atoms, pspots::Array{PsPot_UPF,1} )

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    # From init_run.f90
    nh = zeros(Int64,Nspecies)
    lmaxkb = 0
    for isp in 1:Nspecies
        nh[isp] = 0
        # FIXME: do not add any beta projector if pseudo in 1/r format
        for nb in 1:pspots[isp].Nproj
            l = pspots[isp].proj_l[nb]
            nh[isp] = nh[isp] + 2*l + 1
            lmaxkb = max(lmaxkb, l)
        end
    end
    nhmax = maximum(nh)
    println("nhmax = ", nhmax)

#    for isp in 1:Nspecies
#        for nb in 1:pspots[isp].Nproj
#            l = pspots[isp].proj_l[nb]
#            for m in 1:(2*l + 1)
#                nhmax = nhmax + 1
#            end
#        end
#    end
#    println("nhmax = ", nhmax)


    indv = zeros(Int64,nhmax,Nspecies)         # indes linking  atomic beta's to beta's in the solid
    nhtol = zeros(Int64,nhmax,Nspecies)  # correspondence n <-> angular momentum l
    nhtolm = zeros(Int64,nhmax,Nspecies) # correspondence n <-> combined lm index for (l,m)

    for isp in 1:Nspecies
        ih = 1
        for nb in 1:pspots[isp].Nproj
            l = pspots[isp].proj_l[nb]
            for m in 1:(2*l + 1)
                nhtol[ih,isp] = l
                nhtolm[ih,isp] = l*l+m
                indv[ih,isp] = nb
                ih = ih + 1
            end
        end
    end

    for isp in 1:Nspecies
        println("\nSpecies: ", isp)
        println("nhtol  = ", nhtol[:,isp])
        println("nhtolm = ", nhtolm[:,isp])
        println("indv   = ", indv[:,isp])
    end


    # correspondence beta indexes ih,jh -> composite index ijh
    ijtoh = zeros(Int64,nhmax,nhmax,Nspecies) 
    for isp in 1:Nspecies
        #
        # ijtoh map augmentation channel indexes ih and jh to composite
        # "triangular" index ijh
        ijtoh[:,:,isp] .= -1
        ijv = 0
        for ih in 1:nh[isp], jh in ih:nh[isp]
            ijv = ijv + 1
            ijtoh[ih,jh,isp] = ijv
            ijtoh[jh,ih,isp] = ijv
        end
    end


    indv_ijkb0 = zeros(Int64,Natoms)     # first beta (index in the solid) for each atom 
    ijkb0 = 0
    for isp in 1:Nspecies
        # ijkb0 points to the last beta "in the solid" for atom ia-1
        # i.e. ijkb0+1,.. ijkb0+nh(ityp(ia)) are the nh beta functions of
        #      atom ia in the global list of beta functions (ijkb0=0 for ia=1)
        for ia in 1:Natoms
            if atm2species[ia] == isp
                indv_ijkb0[ia] = ijkb0
                ijkb0 = ijkb0 + nh[isp]
            end
        end
    end
    println("indv_ijkb0 = ", indv_ijkb0)


    # the D functions of the solid
    dvan = zeros(Float64,nhmax,nhmax,Nspecies) 
    for isp in 1:Nspecies
        # From now on the only difference between KB and US pseudopotentials
        # is in the presence of the qq and Q functions.
        #
        # Here we initialize the D of the solid
        for ih in 1:nh[isp], jh in 1:nh[isp]
            cond1 = nhtol[ih,isp] == nhtol[jh,isp]
            cond2 = nhtolm[ih,isp] == nhtolm[jh,isp]
            if cond1 & cond2
                ir = indv[ih,isp]
                is = indv[jh,isp]
                dvan[ih,jh,isp] = pspots[isp].Dion[ir,is]
            end
        end
        display(dvan[:,:,isp]); println()
    end

end