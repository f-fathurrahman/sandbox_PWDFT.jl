using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "GBRV_LDA")


function main_v1()
    atoms = Atoms(xyz_string_frac="""
    2
    
    Si  0.0  0.0  0.0
    Si  0.25  0.25  0.25
    """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631));

    pspots = [ PsPot_UPF(joinpath(DIR_PSP, "si_lda_v1.uspp.F.UPF")) ];
    
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    dual = ecutrho/ecutwfc

    hamOptions = HamiltonianOptions()
    hamOptions.dual = dual
    hamOptions.meshk = [3,3,3]
    Ham = Hamiltonian( atoms, pspots, ecutwfc, hamOptions );
    return
end


function main_v2()
    atoms = Atoms(xyz_string_frac="""
    2
    
    Si  0.0  0.0  0.0
    Si  0.25  0.25  0.25
    """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631));

    pspfiles = [joinpath(DIR_PSP, "si_lda_v1.uspp.F.UPF")];
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, dual=dual, meshk=[3,3,3] );
    return 
end

main_v1()
main_v2()
