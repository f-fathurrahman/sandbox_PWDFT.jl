function gen_eigensystem(
    ispin::Int64, ik::Int64,
    atoms, pw, mt_vars, apwlo_vars, apwlo_ints,
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

    #=
  DO is = 1,nspecies
    ! nuclear charge
    chgzn = chgzn + spzn(is)*natoms(is)
    ! find the maximum number of atomic states
    nstspmax = max(nstspmax,nstsp(is))
    ! compute the electronic charge for each species, as well as the total core and
    ! valence charge
    spze(is) = 0.d0
    chgcr(is) = 0.d0
    DO ist = 1,nstsp(is)
      spze(is) = spze(is) + occsp(ist,is)
      IF( spcore(ist,is) ) THEN 
        chgcr(is) = chgcr(is) + occsp(ist,is)
        nstcr = nstcr + 2*ksp(ist,is)*natoms(is)
      ELSE 
        chgval = chgval + occsp(ist,is)*natoms(is)
      ENDIF 
    ENDDO 
    chgcrtot=chgcrtot + chgcr(is)*natoms(is)
  ENDDO 

    nstfv = int(chgval/2.d0) + nempty + 1
    =#


    # calc eigenvalues and eigenvectors
    evals, evecs = eigen(Hermitian(H), Hermitian(O))

    serialize("H_ispin_ $(ispin)_ik_$(ik).dat", H)
    serialize("O_ispin_$(ispin)_ik_$(ik).dat", O)
    serialize("evals_ispin_$(ispin)_ik_$(ik).dat", evals)
    serialize("evecs_ispin_$(ispin)_ik_$(ik).dat", evecs)

    return
end