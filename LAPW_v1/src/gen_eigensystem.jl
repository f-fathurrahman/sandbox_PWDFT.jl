function gen_eigensystem(
    ispin::Int64, ik::Int64,
    atoms::Atoms, pw::PWGrid, mt_vars::MuffinTins, apwlo_vars::APWLOVars,
    apwlo_ints, elec_chgst,
    nmat, cfunig, vsig
)
    haa = apwlo_ints.haa
    hloa = apwlo_ints.hloa
    hlolo = apwlo_ints.hlolo
    oalo = apwlo_ints.oalo
    ololo = apwlo_ints.ololo

    apwalm = Vector{Array{ComplexF64,3}}(undef, atoms.Natoms)
    for ia in 1:atoms.Natoms
        isp = atoms.atm2species[ia]
        apwordmax = maximum(apwlo_vars.apword[isp])
        Ngk = pw.gvecw.Ngw[ik]
        lmmaxapw = mt_vars.lmmaxapw
        apwalm[ia] = zeros(ComplexF64, Ngk, apwordmax, lmmaxapw)
    end

    calc_match_coeffs!(ik, atoms, pw, mt_vars, apwlo_vars, apwalm)

    # For testing Hamiltonian construction
    #
    H = zeros(ComplexF64, nmat[ik], nmat[ik])

    for ia in 1:atoms.Natoms
        hmlaa!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, haa, H)
    end

    hmlistl!(ik, pw, cfunig, vsig, H)

    for ia in 1:atoms.Natoms
        hmlalo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, hloa, H)
        hmllolo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, hlolo, H)
        # no need to pass apwalm for hmllolo
    end

    O = zeros(ComplexF64, nmat[ik], nmat[ik])
    for ia in 1:atoms.Natoms
        olpaa!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, O)
    end
    olpistl!(ik, pw, cfunig, O)

    for ia in 1:atoms.Natoms
        olpalo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, oalo, O)
        olplolo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, ololo, O)
    end

    # calc eigenvalues and eigenvectors
    evals, evecs = eigen(Hermitian(H), Hermitian(O))

    # FIXME: save full size (?)
    serialize("H_ispin_$(ispin)_ik_$(ik).dat", H)
    serialize("O_ispin_$(ispin)_ik_$(ik).dat", O)
    serialize("evals_ispin_$(ispin)_ik_$(ik).dat", evals)
    serialize("evecs_ispin_$(ispin)_ik_$(ik).dat", evecs)

    
    # Set eigenvalues
    # FIXME: only for the case of nstfv == nstsv
    #@assert elec_chgst.nspinor == 1
    nstsv = elec_chgst.nstsv
    @views elec_chgst.evalsv[1:nstsv,ik] = evals[1:nstsv]

    return
end