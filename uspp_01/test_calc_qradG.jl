using Printf
using SpecialFunctions: sphericalbesselj
using PWDFT

include("loop_calc_qradG.jl")
include("calc_qradG.jl")

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

#
# Main program here
#

atoms = create_atoms_N2H4()
println(atoms)

pw = PWGrid(20.0, atoms.LatVecs, dual=5.0)
println(pw)
    
Nspecies = atoms.Nspecies
#pspfiles = [
#    "/home/efefer/pseudo/GBRV_LDA/n_lda_v1.2.uspp.F.UPF2",
#    "/home/efefer/pseudo/GBRV_LDA/h_lda_v1.4.uspp.F.UPF2"
#]
    
pspfiles = [
    "/home/efefer/pseudo/PSLIB/N.pbe-n-rrkjus_psl.0.1.UPF",
    "/home/efefer/pseudo/PSLIB/H.pbe-rrkjus_psl.0.1.UPF"
]

pspots = Vector{PsPot_UPF}(undef,Nspecies)
for isp in 1:Nspecies
    pspots[isp] = PsPot_UPF( pspfiles[isp] )
    PWDFT._build_prj_interp_table!( pspots[isp], pw )
end

loop_calc_qradG(atoms, pw, pspots)
#qradG = calc_qradG(atoms, pw, pspots)
