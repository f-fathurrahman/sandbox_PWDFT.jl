# Experimenting with USPP

function PsPotNL_UPF(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF}
)

    # Those are parameters (HARDCODED)
    lmaxx  = 3         # max non local angular momentum (l=0 to lmaxx)
    lqmax = 2*lmaxx + 1  # max number of angular momenta of Q

    # calculate the number of beta functions for each atomic type
    nh = zeros(Int64, Nspecies)
    lmaxkb = -1
    for isp in 1:Nspecies
        psp = pspots[isp]
        for nb in 1:psp.Nproj
            nh[isp] = nh[isp] + 2 * psp.proj_l[nb] + 1
            lmaxkb = max(lmaxkb, psp.proj_l[isp])
        end
    end
    println("nh = ", nh)
    println("lmaxkb = ", lmaxkb)


    ap, lpx, lpl = calc_clebsch_gordan(lmaxkb + 1)


    # calculate the maximum number of beta functions
    nhm = maximum(nh)
    # This can be made into a jagged array (vector of vector)
    indv = zeros(Int64, nhm, Nspecies)
    nhtol = zeros(Int64, nhm, Nspecies)
    nhtolm = zeros(Int64, nhm, Nspecies)

    for isp in 1:Nspecies
        ih = 1
        psp = pspots[isp]
        for nb in 1:psp.Nproj
            l = psp.proj_l[nb]
            for m in 1:(2*l+1)
                nhtol[ih,isp] = l
                nhtolm[ih,isp] = l*l + m
                indv[ih,isp] = nb
                ih = ih + 1
            end
        end
    end

    println("indv = ")
    for isp in 1:Nspecies
        println(indv[1:nh[isp],isp])
    end

    println("nhtol = ")
    for isp in 1:Nspecies
        println(nhtol[1:nh[isp],isp])
    end

    println("nhtolm = ")
    for isp in 1:Nspecies
        println(nhtolm[1:nh[isp],isp])
    end

    println("Pass here")

    return

end