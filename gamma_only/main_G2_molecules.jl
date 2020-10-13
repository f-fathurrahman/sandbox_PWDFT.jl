using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("INC_gamma_only.jl")

function do_calc(molname; gamma_only=true)

    Random.seed!(1234)

    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", molname*".xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    
    ecutwfc = 20.0
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc, use_xc_internal=true)
    psis = randn_BlochWavefuncGamma(Ham)
    
    if gamma_only
        #@time KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=200 )
        @time KS_solve_SCF_potmix!( Ham, psis, startingrhoe=:random,
            NiterMax=200, mix_method="broyden" )
    else
        Ham_ = Hamiltonian(atoms, pspfiles, ecutwfc, use_symmetry=false, use_xc_internal=true)
        psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )
        @time KS_solve_Emin_PCG_dot!( Ham_, psiks, startingrhoe=:random,
            skip_initial_diag=true, NiterMax=200 )
    end

end


function main()
    Nargs = length(ARGS)
    if Nargs >= 1
        molname = ARGS[1]
    else
        molname = "H2O"
    end
    for i in 1:1
        do_calc(molname, gamma_only=true)
        #do_calc(molname, gamma_only=false)
    end
end

main()
