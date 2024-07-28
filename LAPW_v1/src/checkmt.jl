function checkmt!(atoms, mt_vars::MuffinTins; epslat=1e-6)
    checkmt!(atoms, mt_vars.rmt, rmtdelta=mt_vars.rmtdelta, epslat=epslat)
    return
end

function checkmt!(atoms::Atoms, specs_info::Vector{SpeciesInfo}; rmtdelta=0.05, epslat=1e-6 )
    Nspecies = atoms.Nspecies
    rmt = zeros(Float64, Nspecies)
    for isp in 1:Nspecies
        rmt[isp] = specs_info[isp].rmt
    end
    checkmt!(atoms, rmt, rmtdelta=rmtdelta, epslat=epslat)
    # Reassign potentially new value
    for isp in 1:Nspecies
        specs_info[isp].rmt = rmt[isp]
    end
    return
end

function checkmt!( atoms::Atoms, rmt::Vector{Float64}; rmtdelta=0.05, epslat=1e-6 )

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols

    rmt0 = zeros(Float64,Nspecies)
    rmt0[1:Nspecies] = rmt[1:Nspecies]

    while true
        # find the minimum distance between muffin-tin surfaces
        dmin, is, js = mtdmin(atoms, rmt)
        # adjust muffin-tin radii if required
        if dmin < (rmtdelta - epslat)
            println("Adjusting rmt")
            t1 = rmt[is] + rmt[js]
            t2 = (t1 + dmin - rmtdelta)/t1
            rmt[is] = rmt[is]*t2
            if is != js
                rmt[js] = rmt[js]*t2
            end
        else
            break
        end
    end

    for is in 1:Nspecies
        if rmt[is] < 0.25
            @printf("Error(checkmt): muffin-tin radius too small for species %d %s\n", is, spsymb[is])
            @printf("Radius : %18.10f\n", rmt[is])
        end    
        if rmt[is] < rmt0[is]
            @printf("Info(checkmt): reduced muffin-tin radius of species %3d %s", is, spsymb[is])
            @printf(" is reduced from %f to %f\n", rmt0[is], rmt[is])
        end
    end

    return
end
