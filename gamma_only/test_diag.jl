using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("INC_gamma_only.jl")
include("diag_Emin_PCG_gamma.jl")

function do_calc(molname; gamma_only=true)

    Random.seed!(1234)

    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", molname*".xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    
    ecutwfc = 15.0
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc)
    psis = randn_BlochWavefuncGamma(Ham)

    Npoints = prod(Ham.pw.Ns)
    Rhoe = zeros(Float64,Npoints,1)

    if gamma_only
        calc_rhoe!(Ham, psis, Rhoe)
        update!(Ham, Rhoe)
        evals = diag_Emin_PCG!(Ham, psis, verbose=true)
    else
        Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
        psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )
        calc_rhoe!(Ham_, psiks, Rhoe)
        update!(Ham_, Rhoe)
        evals = diag_Emin_PCG!(Ham_, psiks, verbose=true)
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
        do_calc(molname, gamma_only=false)
    end
end

main()
