function gensfacgp!(atoms, Ngp, vgpc, sfacgp)
    # XXX Ngp obtained from size(vgpc,2) ?
    for ia in 1:Natoms
        v1 = atposc(1,ia,is)
        v2 = atposc(2,ia,is)
        v3=atposc(3,ia,is)
        for igp in 1:Ngp
            t1 = vgpc[1,igp]*v1 + vgpc[2,igp]*v2 + vgpc[3,igp]*v3
            sfacgp[igp,ia] = cos(t1) + im*sin(t1)
        end
    end
    return
end