using PWDFT

include("calc_rhoe_core.jl")

function main()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = ["/home/efefer/pseudo/ONCV_v0.4.1/nc-sr-04_pw_standard/Si.upf"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    # This should be checked outside this function
    need_rhoe_core = false
    for psp in Ham.pspots
        if psp.is_nlcc
            need_rhoe_core = true
        end
    end

    if need_rhoe_core
        rhoe_core = calc_rhoe_core(Ham.atoms, Ham.pw, Ham.pspots)
    end

    return
end

main()