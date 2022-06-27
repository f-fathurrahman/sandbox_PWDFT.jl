# Experimenting with USPP

struct PsPotNL_UPF
    lmaxx::Int64
    lqmax::Int64
    lmaxkb::Int64
    nh::Vector{Int64}
    nhm::Int64
    ap::Array{Float64,3}
    lpx::Array{Int64,2}
    lpl::Array{Int64,3}
    indv::Array{Int64,2}
    nhtol::Array{Int64,2}
    nhtolm::Array{Int64,2}
    indv_ijkb0::Vector{Int64}
end


function PsPotNL_UPF(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF}
)

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

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
            lmaxkb = max(lmaxkb, psp.proj_l[nb])
        end
    end
    println("nh = ", nh)
    println("lmaxkb = ", lmaxkb)


    ap, lpx, lpl = calc_clebsch_gordan(lmaxkb + 1)


    # calculate the maximum number of beta functions
    nhm = maximum(nh)
    # Some helper indices, for mapping various stuffs
    # Some of these can be made into a jagged array (vector of vector)
    indv = zeros(Int64, nhm, Nspecies)
    nhtol = zeros(Int64, nhm, Nspecies)
    nhtolm = zeros(Int64, nhm, Nspecies)
    indv_ijkb0 = zeros(Int64, Natoms)

    ijkb0 = 0
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
        # ijkb0 points to the last beta "in the solid" for atom ia-1
        # i.e. ijkb0+1,.. ijkb0+nh(ityp(ia)) are the nh beta functions of
        #      atom ia in the global list of beta functions (ijkb0=0 for ia=1)
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            indv_ijkb0[ia] = ijkb0
            ijkb0 = ijkb0 + nh[isp]
        end
    end
    # TODO: Extract lm -> (l,m)

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

    println("nhtolm = ")
    for ia in 1:Natoms
        println(ia, " ", indv_ijkb0[ia])
    end

    return PsPotNL_UPF(
        lmaxx, lqmax, lmaxkb,
        nh, nhm, ap, lpx, lpl,
        indv, nhtol, nhtolm, indv_ijkb0
    )

end



import Base: show
function show( io::IO, psp::PsPotNL_UPF )
    println("PsPotNL_UPF:")
    println("lmaxx = ", psp.lmaxx)
    println("lqmax = ", psp.lqmax)
    println("lmaxkb = ", psp.lmaxkb)
    return
end