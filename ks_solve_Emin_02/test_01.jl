using PWDFT
using Random
using LinearAlgebra
using Printf

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials")

include("KS_solve_Emin_PCG_01.jl")

function create_atoms_N2H4()
    atoms = Atoms(xyz_string="""
6

N       5.94821400       6.81171100       5.22639100
N       5.94821400       5.37379300       5.22639100
H       6.15929600       7.18550400       6.15196500
H       5.00000000       7.09777800       5.00000000
H       5.73713200       5.00000000       6.15196500
H       6.89642800       5.08772600       5.00000000
""", LatVecs=gen_lattice_sc(16.0))
    return atoms
end

function create_Ham_N2H4_oncv()
    atoms = create_atoms_N2H4()
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "N.upf")),
        PsPot_UPF(joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "H.upf"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    options = HamiltonianOptions()
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end



function create_Ham_N2H4_gbrv()
    atoms = create_atoms_N2H4()
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP, "GBRV_LDA", "n_lda_v1.2.uspp.F.UPF")),
        PsPot_UPF(joinpath(DIR_PSP, "GBRV_LDA", "h_lda_v1.4.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end

function create_Ham_N2H4_paw_jth()
    atoms = create_atoms_N2H4()
    pspots = [
        PsPot_UPF("/home/efefer/pseudo/PAW_JTH_LDA/N.upf"),
        PsPot_UPF("/home/efefer/pseudo/PAW_JTH_LDA/H.upf")
    ]
    ecutwfc = 15.0
    ecutrho = 60.0
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end


function main()
    #Ham = create_Ham_N2H4_gbrv()
    Ham = create_Ham_N2H4_paw_jth()
    psiks = rand_BlochWavefunc(Ham)
    KS_solve_Emin_PCG_01!(Ham, psiks)
    #electrons_scf!(Ham, psiks)
end

main()
