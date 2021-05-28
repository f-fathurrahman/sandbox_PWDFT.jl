function potcoul!(
    atoms, atsp_vars,
    mt_vars, pw, 
    rhomt, rhoir
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
        ssz = ssz + sum(rhomt[ia])
    end

    println(rhomt[1][5])
    println(zrhomt[1][5])

    println("sum(rhomt)  = ", ss)
    println("sum(zrhomt) = ", ssz)

    return
end