function potcoul!(
    atoms, atsp_vars,
    mt_vars, pw, 
    rhomt, rhoir,
    vclmt, vclir
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    npmt = mt_vars.npmt
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti

    # Complex rhomt
    zrhomt = Vector{Vector{ComplexF64}}(undef,Natoms)
    ss = 0.0
    ssz = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        zrhomt[ia] = zeros(ComplexF64, npmt[isp])
        r_to_zf_mt!( mt_vars, nrmt[isp], nrmti[isp], rhomt[ia], zrhomt[ia] )
        ss = ss + sum(rhomt[ia])
        ssz = ssz + sum(zrhomt[ia])
    end
    #println("After r_to_zf_mt:")
    #println("sum(rhomt)  = ", ss)
    #println("sum(zrhomt) = ", ssz)

    # FIXME
    zvclmt = Vector{Vector{ComplexF64}}(undef,Natoms)
    ss = 0.0 + im*0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        zvclmt[ia] = zeros(ComplexF64, npmt[isp])
        #
        #println("sum zrhomt[ia] before = ", sum(zrhomt[ia]))
        #println("sum zvclmt[ia] before = ", sum(zvclmt[ia]))
        #
        zpotclmt!( mt_vars, isp, zrhomt[ia], zvclmt[ia] )
        #
        #println("sum zrhomt[ia] after = ", sum(zrhomt[ia]))
        #println("sum zvclmt[ia] after = ", sum(zvclmt[ia]))
        ss += sum(zvclmt[ia])
    end
    #println("\nAfter zpotclmt:")
    #println("sum(zvclmt)  = ", ss)


    # add the nuclear monopole potentials
    vcln = atsp_vars.vcln
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    #println("typeof vcln: ", typeof(vcln))
    #println("size vcln = ", size(vcln))
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        i = 1
        for ir in 1:nri
            zvclmt[ia][i] += vcln[isp][ir]
            i = i + lmmaxi
        end
        for ir in nri+1:nr
            zvclmt[ia][i] += vcln[isp][ir]
            i = i + lmmaxo
        end
    end

    ss = 0.0 + im*0.0
    for ia in 1:Natoms
        ss = ss + sum(zvclmt[ia])
    end
    #println("\nAfter adding vcln:")
    #println("sum(zvclmt)  = ", ss)
    #println(sum(sum.(zvclmt)))

    # store real interstitial charge density in complex array
    Npoints = prod(pw.Ns)
    zrhoir = zeros(ComplexF64, Npoints)
    zrhoir[:] .= rhoir[:]

    # solve Poisson's equation in the entire unit cell
    zvclir = zeros(ComplexF64, Npoints)
    zpotcoul!( atoms, atsp_vars, mt_vars, pw, zrhoir, zvclmt, zvclir )


    # convert complex muffin-tin potential to real spherical harmonic expansion
    for ia in 1:Natoms
        isp = atm2species[ia]
        z_to_rf_mt!( mt_vars, nrmt[isp], nrmti[isp], zvclmt[ia], vclmt[ia])
    end
  
    # store complex interstitial potential in real array
    @views vclir[:] = real(zvclir[:])

    #println("sum vclmt = ", sum(sum.(vclmt)))
    #println("sum vclir = ", sum(abs.(vclir)))

    return
end