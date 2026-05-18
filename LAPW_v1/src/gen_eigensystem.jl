function gen_eigensystem!(
    ik::Int64,
    atoms::Atoms, atsp_vars, pw::PWGrid, mt_vars::MuffinTins, apwlo_vars::APWLOVars,
    apwlo_ints, elec_chgst,
    nmat, cfunig, vsig;
    bsmt=nothing, bsir=nothing, ndmag=0,
    tmpdir = "./tmp"
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
    calc_match_coeffs!(ik, atoms, atsp_vars, pw, mt_vars, apwlo_vars, apwalm)
    #for ia in 1:atoms.Natoms
    #    println("gen_eigensystem: ik = $ik, ia = $ia, sum abs apwalm = ", sum(abs.(apwalm[ia])))
    #end

    # For testing Hamiltonian construction
    H = zeros(ComplexF64, nmat[ik], nmat[ik])

    for ia in 1:atoms.Natoms
        hmlaa!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, haa, H)
        #println("hmlaa, ia=$ia sum abs H = ", sum(abs.(H)))
    end
    #println("Before hmlist, sum abs H = ", sum(abs.(H)))
    hmlistl!(ik, pw, cfunig, vsig, H)
    #println("After hmlist, sum abs H = ", sum(abs.(H)))

    for ia in 1:atoms.Natoms
        hmlalo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, hloa, H)
        hmllolo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, hlolo, H)
        # no need to pass apwalm for hmllolo
    end

    O = zeros(ComplexF64, nmat[ik], nmat[ik])
    for ia in 1:atoms.Natoms
        olpaa!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, O)
        #println("olpaa, ia=$ia sum abs O = ", sum(abs.(O)))
    end
    #println("Before olpist, sum abs O = ", sum(abs.(O)))
    olpistl!(ik, pw, cfunig, O)
    #println("After olpist, sum abs O = ", sum(abs.(O)))

    for ia in 1:atoms.Natoms
        olpalo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, oalo, O)
        olplolo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, ololo, O)
    end

    #println("sum abs H = ", sum(abs.(H)))
    #println("sum abs O = ", sum(abs.(O)))
    # calc eigenvalues and eigenvectors
    evals, evecs = eigen(Hermitian(H), Hermitian(O))

    nstsv = elec_chgst.nstsv
    nstfv = elec_chgst.nstfv

    # FIXME: save full size (?)
    #serialize("H_ik_$(ik).dat", H)
    #serialize("O_ik_$(ik).dat", O)
    serialize(joinpath(tmpdir, "evals_1st_ik_$(ik).dat"), evals[1:nstfv])
    serialize(joinpath(tmpdir, "evecs_1st_ik_$(ik).dat"), evecs[:,1:nstfv])
    
    # Set eigenvalues
    # FIXME: only for the case of nstfv == nstsv
    #@assert elec_chgst.nspinor == 1
    #@views elec_chgst.evalsv[1:nstsv,ik] = evals[1:nstsv]

    gen_eigensystem_2nd!(
        ik,
        atoms, pw, mt_vars, apwlo_vars,
        apwalm, elec_chgst, bsmt, bsir, ndmag,
        tmpdir = tmpdir
    )

    return
end