using LinearAlgebra
using Printf
using Random

import PWDFT
using PWDFT: Hamiltonian, PsPot_GTH, Atoms, PWGrid, ANG2BOHR
using PWDFT: gen_lattice_fcc, gen_lattice_sc
using PWDFT: KS_solve_SCF!, R_to_G, eval_Vloc_G, calc_strfact, G_to_R, KS_solve_Emin_PCG!

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")
include("calc_forces_Ps_loc.jl")

function main()

    #Ham = create_Ham_H2()
    #Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_Si_fcc(xcfunc="PBE")
    Ham = create_Ham_GaAs_v1()
    #Ham = create_Ham_Si8()
    #Ham = create_Ham_GaAs_v2()
    #Ham = create_Ham_CO()

    Random.seed!(1234)
    KS_solve_Emin_PCG!(Ham, etot_conv_thr=1e-8)

    Natoms = Ham.atoms.Natoms
    atsymbs = Ham.atoms.atsymbs

    println("F_Ps_loc: ")
    F_Ps_loc = calc_forces_Ps_loc(Ham)*2.0
    @time F_Ps_loc = calc_forces_Ps_loc(Ham)*2.0
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

    println("PWDFT F_Ps_loc: ")
    F_Ps_loc = PWDFT.calc_forces_Ps_loc(Ham)*2.0
    @time F_Ps_loc = PWDFT.calc_forces_Ps_loc(Ham)*2.0
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

    println("")
    println("Using finite difference")
    F_Ps_loc = calc_forces_Ps_loc_finite_diff( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
    @time F_Ps_loc = calc_forces_Ps_loc_finite_diff( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

end


main()