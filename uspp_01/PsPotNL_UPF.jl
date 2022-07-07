# Experimenting with USPP

mutable struct PsPotNL_UPF
    lmaxx::Int64
    lqmax::Int64
    lmaxkb::Int64
    nh::Vector{Int64}
    nhm::Int64
    nkb::Int64
    ap::Array{Float64,3}
    lpx::Array{Int64,2}
    lpl::Array{Int64,3}
    indv::Array{Int64,2}
    nhtol::Array{Int64,2}
    nhtolm::Array{Int64,2}
    indv_ijkb0::Vector{Int64}
    Dvan::Array{Float64,3}
    qradG::Union{Vector{Array{Float64,3}},Nothing}
    qq_nt::Union{Array{Float64,3},Nothing}
    qq_at::Union{Array{Float64,3},Nothing}
    betaNL::Vector{Matrix{ComplexF64}}
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

    nkb = 0
    for ia in 1:Natoms
       isp = atm2species[ia]
       nkb = nkb + nh[isp]
    end

    ap, lpx, lpl = calc_clebsch_gordan(lmaxkb + 1)


    # calculate the maximum number of beta functions
    nhm = maximum(nh)
    # Some helper indices, for mapping various stuffs
    # Some of these can be made into a jagged array (vector of vector)
    indv = zeros(Int64, nhm, Nspecies)
    nhtol = zeros(Int64, nhm, Nspecies)
    nhtolm = zeros(Int64, nhm, Nspecies)
    indv_ijkb0 = zeros(Int64, Natoms)
    Dvan = zeros(Float64, nhm, nhm, Nspecies)

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

        for ih in 1:nh[isp], jh in 1:nh[isp]
            cond_l  = nhtol[ih,isp] == nhtol[jh,isp]
            cond_lm = nhtolm[ih,isp] == nhtolm[jh,isp]
            if cond_l && cond_lm
                ihs = indv[ih,isp]
                jhs = indv[jh,isp]
                Dvan[ih,jh,isp] = psp.Dion[ihs,jhs]
            end
        end
    end
    # TODO: Extract lm -> (l,m)

    tmp_uspp = zeros(Bool,Nspecies)
    for isp in 1:Nspecies
        tmp_uspp[isp] = pspots[isp].is_ultrasoft
    end
    #
    if all(tmp_uspp)
        qradG, qq_nt, qq_at = _prepare_aug_charges(
            atoms, pw, pspots, lmaxkb, nhm, nh, indv, nhtolm, lpl, lpx, ap
        )
    else
        qradG = nothing
        qq_nt = nothing
        qq_at = nothing
    end

    Nkpt = pw.gvecw.kpoints.Nkpt
    betaNL = Vector{Matrix{ComplexF64}}(undef,Nkpt)
    for ik in 1:Nkpt
        betaNL[ik] = zeros(ComplexF64, pw.gvecw.Ngw[ik], nkb)
        _init_Vnl_KB!(
            ik, atoms, pw, pspots,
            lmaxkb, nh, nhm, nhtol, nhtolm, indv,
            betaNL[ik]
        )
    end

    return PsPotNL_UPF(
        lmaxx, lqmax, lmaxkb,
        nh, nhm, nkb, ap, lpx, lpl,
        indv, nhtol, nhtolm, indv_ijkb0,
        Dvan, qradG, qq_nt, qq_at,
        betaNL
    )

end


function _prepare_aug_charges(
    atoms, pw, pspots,
    lmaxkb, nhm, nh, indv, nhtolm, lpl, lpx, ap
)
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    qradG = calc_qradG(pw, pspots)

    # Compute the qq coefficients by integrating the Q.
    # The qq are the g=0 components of Q

    # FIXME: use Vector{Matrix} instead of 3d array
    qq_nt = zeros(Float64, nhm, nhm, Nspecies) # qq_nt, need qvan2
    qq_at = zeros(Float64, nhm, nhm, Natoms)

    G0 = zeros(3,1)
    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, 1, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    Ylm_real_qe!(_lmax, G0, ylmk0) # Ylm_real_qe accept l value starting from 0
    qgm = zeros(ComplexF64, 1)
    for isp in 1:Nspecies
        for ih in 1:nh[isp], jh in ih:nh[isp]
            qvan2!( indv, nhtolm, lpl, lpx, ap, qradG, ih, jh, isp, [0.0], ylmk0, qgm )
            qq_nt[ih,jh,isp] = pw.CellVolume * real(qgm[1])
            qq_nt[jh,ih,isp] = pw.CellVolume * real(qgm[1])
        end
    end

    # finally we set the atomic specific qq_at matrices
    for ia in 1:Natoms
        qq_at[:,:,ia] = qq_nt[:,:,atm2species[ia]]
    end

    return qradG, qq_nt, qq_at
end

include("init_Vnl_KB.jl")


import Base: show
function show( io::IO, pspotNL::PsPotNL_UPF )
    
    println("------------")
    println("PsPotNL_UPF:")
    println("------------")
    
    println("lmaxx  = ", pspotNL.lmaxx)
    println("lqmax  = ", pspotNL.lqmax)
    println("lmaxkb = ", pspotNL.lmaxkb)
    println("nkb    = ", pspotNL.nkb)

    nh = pspotNL.nh
    indv = pspotNL.indv
    nhtol = pspotNL.nhtol
    nhtolm = pspotNL.nhtolm
    indv_ijkb0 = pspotNL.indv_ijkb0
    Dvan = pspotNL.Dvan
    qq_nt = pspotNL.qq_nt

    Nspecies = size(indv,2)
    Natoms = size(indv_ijkb0,1)

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

    println("Dvan = ")
    for isp in 1:Nspecies
        println("Dvan isp = ", isp)
        display(Dvan[:,:,isp]); println()
    end


    for isp in 1:Nspecies
        if qq_nt != nothing
            println("qq_nt")
            display(qq_nt[1:nh[isp],1:nh[isp],isp]); println();
        end
    end

    return
end