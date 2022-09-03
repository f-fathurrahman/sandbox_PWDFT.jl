using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("create_Ham.jl")

function test_main()

    #Ham = create_Ham_H2()
    #Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_Si_fcc(xcfunc="PBE")
    Ham = create_Ham_GaAs_v1()
    #Ham = create_Ham_Si8()
    #Ham = create_Ham_GaAs_v2()
    #Ham = create_Ham_CO()
    #println(Ham)

    Random.seed!(1234)

    psiks = rand_BlochWavefunc(Ham)
    #KS_solve_SCF!( Ham, psiks, mix_method="anderson", etot_conv_thr=1e-6 )
    #KS_solve_SCF!( Ham, psiks, mix_method="broyden", etot_conv_thr=1e-8 )
    KS_solve_Emin_PCG!( Ham, psiks, etot_conv_thr=1e-8 )

    atoms = Ham.atoms
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    F_NN = calc_forces_NN( Ham.pw, Ham.atoms )*2.0
    println("NN forces: (in Ry/bohr)")
    @time F_NN = calc_forces_NN( Ham.pw, Ham.atoms )*2.0
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

    F_Ps_loc = calc_forces_Ps_loc( Ham )*2.0
    println("Ps loc forces: (in Ry/bohr)")
    @time F_Ps_loc = calc_forces_Ps_loc( Ham )*2.0
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

    F_Ps_nloc = calc_forces_Ps_nloc( Ham, psiks )*2.0
    println("Ps nloc forces: (in Ry/bohr)")
    @time F_Ps_nloc = calc_forces_Ps_nloc( Ham, psiks )*2.0
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end

    F_total = F_NN + F_Ps_loc + F_Ps_nloc
    println("\nTotal forces: (in Ry unit)")
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_total[1,ia], F_total[2,ia], F_total[3,ia] )
    end

end


test_main()
