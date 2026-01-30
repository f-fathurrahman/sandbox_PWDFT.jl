
function calc_engynn(atoms, atsp_vars, pw, mt_vars)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    vcln = atsp_vars.vcln

    npmt = mt_vars.npmt
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo

    Npoints = prod(pw.Ns)

    zvclmt = Vector{Vector{ComplexF64}}(undef, Natoms)

    # generate the nuclear monopole potentials
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        zvclmt[ia] = zeros(ComplexF64, npmt[isp])
        i = 1
        for ir in 1:nri
            zvclmt[ia][i] = vcln[isp][ir]
            i += lmmaxi
        end
        for ir in (nri+1):nr
            zvclmt[ia][i] = vcln[isp][ir]
            i += lmmaxo
        end
    end
    zrhoir = zeros(ComplexF64, Npoints)
    zvclir = zeros(ComplexF64, Npoints)
    # solve the complex Poisson's equation
    zpotcoul!( atoms, atsp_vars, mt_vars, pw, zrhoir, zvclmt, zvclir )
    
    # compute the nuclear-nuclear energy
    E_nn = 0.0
    y00 = 0.5/sqrt(pi)
    Zatoms = -PWDFT.get_Zatoms(atoms) # use minus sign
    for ia in 1:Natoms
        isp = atm2species[ia]
        t1 = ( real(zvclmt[ia][1]) - vcln[isp][1] )*y00
        E_nn += Zatoms[isp]*t1
    end
    E_nn *= 0.5
    return E_nn
end