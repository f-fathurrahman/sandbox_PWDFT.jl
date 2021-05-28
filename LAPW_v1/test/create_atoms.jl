function create_H2O()
    LatVecs = zeros(3,3)
    A = 6.0
    LatVecs[1,:] = [A, 0.0, 0.0]
    LatVecs[2,:] = [0.0, A, 0.0]
    LatVecs[3,:] = [0.0, 0.0, A]

    atoms = Atoms(xyz_string_frac="""
    3

    O   0.0000000000      0.0000000000      0.0000000000
    H   0.3018333333      0.0000000000      0.0000000000
    H   0.9246820068      0.2922850682      0.0000000000
    """, in_bohr=true, LatVecs=LatVecs)

    return atoms
end

function create_Si_fcc()
    LatVecs = zeros(3,3)
    A = 5.13
    LatVecs[1,:] = [A, A, 0.0]
    LatVecs[2,:] = [A, 0.0, A]
    LatVecs[3,:] = [0.0, A, A]

    atoms = Atoms(xyz_string_frac="""
    2

    Si  0.0  0.0  0.0
    Si  0.25 0.25 0.25
    """, in_bohr=true, LatVecs=LatVecs)

    return atoms
end

function create_SiPt_fcc()
    LatVecs = zeros(3,3)
    A = 5.13
    LatVecs[1,:] = [A, A, 0.0]
    LatVecs[2,:] = [A, 0.0, A]
    LatVecs[3,:] = [0.0, A, A]

    atoms = Atoms(xyz_string_frac="""
    2

    Si  0.0  0.0  0.0
    Pt  0.25 0.25 0.25
    """, in_bohr=true, LatVecs=LatVecs)

    return atoms
end

function create_Si_atom()
    LatVecs = zeros(3,3)
    A = 6.0
    LatVecs[1,:] = [A, 0.0, 0.0]
    LatVecs[2,:] = [0.0, A, 0.0]
    LatVecs[3,:] = [0.0, 0.0, A]

    atoms = Atoms(xyz_string_frac="""
    1

    Si  0.0  0.0  0.0
    """, in_bohr=true, LatVecs=LatVecs)

    return atoms
end