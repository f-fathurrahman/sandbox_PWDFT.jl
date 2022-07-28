function init_Ham_Si_dimer()
    # Atoms
    atoms = Atoms(xyz_string="""
    2

    Si  3.505  2.5   2.5
    Si  0.0    2.5   2.5
    """, in_bohr=true, LatVecs = gen_lattice_sc(10.0) )

    write_xsf("TEMP_Si_dimer.xsf", atoms)

    pspfiles = ["Si-q4_mod.gth"] # use a modified pspot
    #pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 10.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    Ham.atoms.masses[:] = [28.085]*AMU_AU
    return Ham
end

function init_Ham_H2O()
    ecutwfc = 15.0
    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", "H2O.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [16.0, 2.0]*AMU_AU
    #
    return Ham
end


function init_Ham_CO2()
    ecutwfc = 15.0
    #filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", "CO2.xyz")
    filename = joinpath("../md_01/CO2.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [14.0, 16.0]*AMU_AU

    return Ham
end


function init_Ham_Si8()
    ecutwfc = 10.0
    atoms = Atoms(xyz_string="""
    8

    Si       0.1000000000       0.0000000000       0.0000000000
    Si       7.6969017521       7.6969017521       2.5656339174
    Si       5.1312678348       0.0000000000       5.1312678348
    Si       7.6969017521       2.5656339174       7.6969017521
    Si       0.0000000000       5.1312678348       5.1312678348
    Si       2.5656339174       2.5656339174       2.5656339174
    Si       2.5656339174       7.6969017521       7.6969017521
    Si       5.1312678348       5.1312678348       0.0000000000
    """, in_bohr=true, LatVecs=gen_lattice_sc(10.2625356695))
    
    write_xsf("TEMP_Si8.xsf", atoms)

    #pspfiles = get_default_psp(atoms)
    pspfiles = ["Si-q4_mod.gth"] # use a modified pspot
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [28.0855]*AMU_AU
    return Ham
end