function _init_Vnl_KB!(
    ik::Int64,
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    pspotNL::PsPotNL_UPF,
    Vnl_KB::Array{ComplexF64,2}
)
    _init_Vnl_KB!(
        ik, atoms, pw, pspots,
        pspotNL.lmaxkb, pspotNL.nh, pspotNL.nhm,
        pspotNL.nhtol, pspotNL.nhtolm, pspotNL.indv,
        Vnl_KB
    )
    return
end



# From init_us_2
function _init_Vnl_KB!(
    ik::Int64,
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    lmaxkb, nh, nhm, nhtol, nhtolm, indv,
    Vnl_KB::Array{ComplexF64,2}
)

    Ngw = pw.gvecw.Ngw
    @assert size(Vnl_KB,1) == Ngw[ik]

    idx_gw2g = pw.gvecw.idx_gw2g
    k = pw.gvecw.kpoints.k
    G = pw.gvec.G

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
                # No need to offset nhtol by 1. It is already physics's l (start from l=0)
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