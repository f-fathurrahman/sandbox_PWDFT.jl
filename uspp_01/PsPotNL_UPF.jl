# Experimenting with USPP

struct PsPotNL_UPF
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
    qradG::Vector{Array{Float64,3}}
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
    qq_nt = zeros(Float64, nhm, nhm, Nspecies) # qq_nt, need qvan2


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

    qradG = calc_qradG(pw, pspots) # FIXME: include in PsPotNL_UPF


#=
    G0 = [0.0, 0.0, 0.0]
    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, 1, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    Ylm_real_qe!(_lmax, G0, ylmk0) # Ylm_real_qe accept l value starting from 0
    qgm = zeros(ComplexF64, 1)
    for isp in 1:Nspecies
        for ih in 1:nh[isp], for jh in ih:nh[isp]
            qvan2!( indv, nhtolm, lpl, lpx, ap, 1, ih, jh, isp, G0, ylmk0, qgm )
            qq_nt[ih,jh,nt] = pw.CellVolume * real(qgm[1])
            qq_nt[jh,ih,nt] = pw.CellVolume * real(qgm[1])
        end
    end
=#


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
        println("Dion: ")
        display(pspots[isp].Dion); println()
    end


    return PsPotNL_UPF(
        lmaxx, lqmax, lmaxkb,
        nh, nhm, nkb, ap, lpx, lpl,
        indv, nhtol, nhtolm, indv_ijkb0,
        Dvan, qradG
    )

end


# From init_us_2
function _init_Vnl_KB!(
    ik::Int64,
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    pspotNL::PsPotNL_UPF,
    Vnl_KB::Array{ComplexF64,2}
)

    Ngw = pw.gvecw.Ngw
    @assert size(Vnl_KB,1) == Ngw[ik]

    idx_gw2g = pw.gvecw.idx_gw2g
    k = pw.gvecw.kpoints.k
    G = pw.gvec.G

    lmaxkb = pspotNL.lmaxkb
    nhtolm = pspotNL.nhtolm
    nhtol = pspotNL.nhtol
    nhm = pspotNL.nhm
    nh = pspotNL.nh
    indv = pspotNL.indv

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atpos = atoms.positions
    atm2species = atoms.atm2species

    Gk = zeros(Float64, 3, Ngw[ik])
    Gk2 = zeros(Float64, Ngw[ik])
    # 
    # igk: index for Gk
    # ig: index for G
    for igk in 1:Ngw[ik]
        ig = idx_gw2g[ik][igk] # index of Gk in G
        Gk[1,igk] = G[1,ig] + k[1,ik]
        Gk[2,igk] = G[2,ig] + k[2,ik]
        Gk[3,igk] = G[3,ig] + k[3,ik]
        Gk2[igk] = Gk[1,igk]^2 +  Gk[2,igk]^2 + Gk[3,igk]^2
        # need Gk2? it is also calculated in Ylm_real_qe!
    end

    ylm = zeros(Float64, Ngw[ik], (lmaxkb+1)^2)
    # Ylm_real_qe accept l value starting from 0 (the actual 'physics' angular momentum number)
    Ylm_real_qe!(lmaxkb, Gk, ylm)


    vq = zeros(Float64, Ngw[ik])
    vkb1 = zeros(Float64, Ngw[ik], nhm)
    Sf = zeros(ComplexF64, Ngw[ik])
    #
    # |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
    jkb = 0 # index of beta functions
    dq = 0.01 # HARDCODED
    #    
    for isp in 1:Nspecies
        #
        psp = pspots[isp]
        tab = psp.prj_interp_table
        println("isp = ", isp)
        println("psp.Nproj = ", psp.Nproj)
        # calculate beta in G-space using an interpolation table:
        # f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
        for ibeta in 1:psp.Nproj
            #
            for igk in 1:Ngw[ik]
                Gm = sqrt(Gk2[igk])
                px = Gm/dq - floor(Int64, Gm/dq )
                ux = 1.0 - px
                vx = 2.0 - px
                wx = 3.0 - px
                i0 = floor(Int64, Gm/dq ) + 1
                i1 = i0 + 1
                i2 = i0 + 2
                i3 = i0 + 3
                vq[igk] = tab[i0,ibeta] * ux * vx * wx / 6.0 +
                          tab[i1,ibeta] * px * vx * wx / 2.0 -
                          tab[i2,ibeta] * px * ux * wx / 2.0 +
                          tab[i3,ibeta] * px * ux * vx / 6.0
            end

      
            # add spherical harmonic part  (Y_lm(q)*f_l(q)) 
            for ih in 1:nh[isp]
                if ibeta == indv[ih,isp]        
                    lm = nhtolm[ih,isp]
                    for igk in 1:Ngw[ik]
                        vkb1[igk,ih] = ylm[igk,lm] * vq[igk]
                    end
                end
            end

        end

        println("sum vkb1 = ", sum(vkb1[:,1:nh[isp]]))

        # vkb1 contains all betas including angular part for type nt
        # now add the structure factor and factor (-i)^l

        # ordering: first all betas for atoms of type 1
        #           then  all betas for atoms of type 2  and so on
        #
        for ia in 1:Natoms
            # skip if this is not the current species index
            if atm2species[ia] != isp
                continue
            end
            #
            for igk in 1:Ngw[ik]
                # XXX use dot product here?
                GkX = atpos[1,ia]*Gk[1,igk] + atpos[2,ia]*Gk[2,igk] + atpos[3,ia]*Gk[3,igk]
                Sf[igk] = cos(GkX) - im*sin(GkX)
            end
            
            for ih in 1:nh[isp]
                # XXXX idx of KB projectors increment here ...
                jkb = jkb + 1
                #
                pref = (-im)^nhtol[ih,isp]
                for igk in 1:Ngw[ik]
                    Vnl_KB[igk,jkb] = vkb1[igk,ih] * Sf[igk] * pref
                end
                # clean up garbage in the last block
                #DO ig = npw_+1, npwx
                #    vkb_(ig, jkb) = (0.0_DP, 0.0_DP)
                #ENDDO
            end
        end
    
    end

    println("sum Vnl_KB = ", sum(Vnl_KB))
    println("shape Vnl_KB = ", size(Vnl_KB))

    return

end



import Base: show
function show( io::IO, pspotNL::PsPotNL_UPF )
    println("PsPotNL_UPF:")
    println("lmaxx  = ", pspotNL.lmaxx)
    println("lqmax  = ", pspotNL.lqmax)
    println("lmaxkb = ", pspotNL.lmaxkb)
    println("nkb    = ", pspotNL.nkb)
    return
end