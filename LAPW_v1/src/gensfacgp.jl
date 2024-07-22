function gensfacgp!(atoms, Ngp, vgpc, sfacgp)
    atposc = atoms.positions
    Natoms = atoms.Natoms
    # XXX Ngp obtained from size(vgpc,2) ?
    for ia in 1:Natoms
        v1 = atposc[1,ia]
        v2 = atposc[2,ia]
        v3 = atposc[3,ia]
        for igp in 1:Ngp
            t1 = vgpc[1,igp]*v1 + vgpc[2,igp]*v2 + vgpc[3,igp]*v3
            sfacgp[igp,ia] = cos(t1) + im*sin(t1)
        end
    end
    return
end