function findsymsite!(sym_vars::SymmetryVars, atoms::Atoms)

    Natoms = atoms.Natoms
    atposc = atoms.positions
    atposl = inv(atoms.LatVecs)*atposc

    nsymsite = zeros(Int64, Natoms)
    lspnsyms = zeros(Int64, 48, Natoms)
    lsplsyms = zeros(Int64, 48, Natoms)
    
    apl = zeros(Float64, 3, Natoms)
    iea = zeros(Int64, Natoms, 48)

    for ia in 1:Natoms
        for ja in 1:Natoms
            apl[:,ja] = atposl[:,ja] - atposl[:,ia]
        end
        @views nsymsite[ia] = findsym!(sym_vars, atoms, apl, apl, lsplsyms[:,ia], lspnsyms[:,ia], iea)
    end

    # update the corresponding variables in sym_vars
    sym_vars.nsymsite = nsymsite
    sym_vars.lspnsyms = lspnsyms
    sym_vars.lsplsyms = lsplsyms

    return

end

