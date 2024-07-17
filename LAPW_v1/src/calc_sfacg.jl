# Related function: gensfacgp
#
# Structure factor, using Elk convention, i.e (+) sign in the exponential
function calc_sfacg(atoms, pw)
    Natoms = atoms.Natoms
    Ng = pw.gvec.Ng
    sfacg = zeros(ComplexF64, Ng, Natoms)
    for ia in 1:Natoms
        v1 = atoms.positions[1,ia]
        v2 = atoms.positions[2,ia]
        v3 = atoms.positions[3,ia]
        for ig in 1:Ng
            Gv1 = pw.gvec.G[1,ig]
            Gv2 = pw.gvec.G[2,ig]
            Gv3 = pw.gvec.G[3,ig]
            x = v1*Gv1 + v2*Gv2 + v3*Gv3
            sfacg[ig,ia] = cos(x) + im*sin(x)
        end
    end 
    return sfacg
end