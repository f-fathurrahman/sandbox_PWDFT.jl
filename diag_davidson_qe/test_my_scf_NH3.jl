using LinearAlgebra
using Random
using Printf
using PWDFT

include("op_S_identity.jl")
include("diag_davidson_qe_v2.jl")
include("my_scf.jl")

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function test_my_scf(β::Float64)
    
    Random.seed!(1234)

    atoms = Atoms(
        xyz_file=joinpath(DIR_PWDFT, "structures", "NH3.xyz"),
        LatVecs=gen_lattice_sc(16.0)
    )
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]

    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, 15.0 )

    psiks = rand_BlochWavefunc(Ham)
    my_scf!(Ham, psiks, betamix=β)

end

test_my_scf(0.2)
@time test_my_scf(0.2)

