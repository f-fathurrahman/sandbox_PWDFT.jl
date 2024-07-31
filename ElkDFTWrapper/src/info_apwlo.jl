# Only print out the important data



function info_apwve()
    #
    nspecies = get_nspecies()
    apword = get_apword() # actual order of APW
    lmaxapw = get_lmaxapw() # actual lmax used in APW expansion
    println()
    println("-----------")
    println("INFO apwve")
    println("-----------")
    #
    apwve = get_apwve()
    # apwve(maxapword,0:maxlapw,maxspecies)
    for isp in 1:nspecies
        println()
        println("Species index = ", isp)
        for l in 0:lmaxapw
            println("APW order for isp=$isp and l=$l = $(apword[l,isp])")
            for io in 1:apword[l,isp]
                println("\tisp=$isp l=$l io=$io apwve=$(apwve[io,l,isp])")
            end
        end
    end
    return
end



function info_apwe0()
    #
    nspecies = get_nspecies()
    apword = get_apword() # actual order of APW
    lmaxapw = get_lmaxapw() # actual lmax used in APW expansion
    println()
    println("----------")
    println("INFO apwe0")
    println("----------")
    #
    apwe0 = get_apwe0()
    # apwe0(maxapword,0:maxlapw,maxspecies)
    for isp in 1:nspecies
        println()
        println("Species index = ", isp)
        for l in 0:lmaxapw
            println("APW order for isp=$isp and l=$l = $(apword[l,isp])")
            for io in 1:apword[l,isp]
                println("\tisp=$isp l=$l io=$io apwe0=$(apwe0[io,l,isp])")
            end
        end
    end
    return
end



function info_apwe()
    #
    nspecies = get_nspecies()
    natoms = get_natoms() # per-species
    natmtot = get_natmtot()
    idxas = get_idxas() # atom within species to global atom index
    apword = get_apword() # actual order of APW
    lmaxapw = get_lmaxapw() # actual lmax used in APW expansion
    println()
    println("---------")
    println("INFO apwe")
    println("---------")
    #
    apwe = get_apwe()
    # apwe(apwordmax, 0:lmaxapw,natmtot)
    for isp in 1:nspecies, ias in 1:natoms[isp]
        ia = idxas[ias,isp]
        println()
        println("Atom index = ", ia)
        for l in 0:lmaxapw
            println("APW order for isp=$isp and l=$l = $(apword[l,isp])")
            for io in 1:apword[l,isp]
                println("\tia=$ia l=$l io=$io apwe=$(apwe[io,l,ia])")
            end
        end
    end
    return
end


function info_lorbve()
    #
    nspecies = get_nspecies()
    lorbord = get_lorbord() # actual order or LO for species and local orbital
    nlorb = get_nlorb() # actual number of LO
    println()
    println("-----------")
    println("INFO lorbve")
    println("-----------")
    #
    lorbve = get_lorbve()
    # lorbve(maxlorbord,maxlorb,maxspecies)
    for isp in 1:nspecies
        println()
        println("Species index = ", isp)
        for ilo in 1:nlorb[isp]
            println("LO order for isp=$isp and ilo=$ilo = $(lorbord[ilo,isp])")
            for io in 1:lorbord[ilo,isp]
                println("\tisp=$isp ilo=$ilo io=$io lorbve=$(lorbve[io,ilo,isp])")
            end
        end
    end
    return
end



function info_lorbe0()
    #
    nspecies = get_nspecies()
    lorbord = get_lorbord() # actual order or LO for species and local orbital
    nlorb = get_nlorb() # actual number of LO
    println()
    println("-----------")
    println("INFO lorbe0")
    println("-----------")
    #
    lorbe0 = get_lorbe0()
    # lorbe0(maxlorbord,maxlorb,maxspecies)
    for isp in 1:nspecies
        println()
        println("Species index = ", isp)
        for ilo in 1:nlorb[isp]
            println("LO order for isp=$isp and ilo=$ilo = $(lorbord[ilo,isp])")
            for io in 1:lorbord[ilo,isp]
                println("\tisp=$isp ilo=$ilo io=$io lorbe0=$(lorbe0[io,ilo,isp])")
            end
        end
    end
    return
end

function info_lorbe()
    #
    nspecies = get_nspecies()
    natoms = get_natoms() # per-species
    natmtot = get_natmtot()
    idxas = get_idxas() # atom within species to global atom index
    lorbord = get_lorbord() # actual order or LO for species and local orbital
    nlorb = get_nlorb() # actual number of LO
    println()
    println("----------")
    println("INFO lorbe")
    println("----------")
    #
    lorbe = get_lorbe()
    # lorbe(lorbordmax,maxlorb,natmtot)
    for isp in 1:nspecies, ias in 1:natoms[isp]
        ia = idxas[ias,isp]
        println()
        println("Atom index = ", ia)
        for ilo in 1:nlorb[isp]
            println("LO order for isp=$isp and ilo=$ilo = $(lorbord[ilo,isp])")
            for io in 1:lorbord[ilo,isp]
                println("\tia=$(ia) ilo=$ilo io=$io lorbe=$(lorbe[io,ilo,ia])")
            end
        end
    end
    return
end
