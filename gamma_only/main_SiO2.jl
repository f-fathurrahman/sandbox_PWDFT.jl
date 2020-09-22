using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("INC_gamma_only.jl")

function do_calc(;gamma_only=true)

    Random.seed!(1234)

    # In Bohr
    A = 9.28990
    B = 9.28990*1.73206
    C = 9.28990*1.09955

    atoms = Atoms(xyz_string="""
    18

     O  3.18829368  14.83237039   1.22882961
     O  7.83231469   6.78704039   1.22882961
     O  2.07443467   5.99537992   4.73758250
     O  6.72031366  14.04231898   4.73758250
     O  3.96307134  11.26989826   7.87860582
     O  8.60802134   3.22295920   7.87860582
     O  3.96307134   4.81915267   9.14625133
     O  8.60802134  12.86448267   9.14625133
     O  3.18736469   1.25668055   5.58029607
     O  7.83324368   9.30201055   5.58029607
     O  2.07536366  10.09206195   2.07358613
     O  6.71938467   2.04673195   2.07358613
    Si  0.28891589   8.04533000   3.40456284
    Si  4.93386589   0.00000000   3.40456284
    Si  2.13389003  12.27717358  -0.04188031
    Si  6.77884003   4.23184358  -0.04188031
    Si  2.13389003   3.81348642   6.85202747
    Si  6.77884003  11.85881642   6.85202747
    """, in_bohr=true, LatVecs=diagm([A,B,C]))
    pspfiles = get_default_psp(atoms)
    
    ecutwfc = 15.0
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc)
    psis = randn_BlochWavefuncGamma(Ham)
    
    if gamma_only
        @time KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=200 )
    else
        Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
        psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )
        @time KS_solve_Emin_PCG_dot!( Ham_, psiks, startingrhoe=:random,
            skip_initial_diag=true, NiterMax=200 )
    end

end


do_calc(gamma_only=true)
