using Printf
using OffsetArrays

using PWDFT

include("create_atoms_N2H4.jl")
include("../ylm_real/Ylm_real_qe.jl")
include("calc_clebsch_gordan.jl")
include("PsPotNL_UPF.jl")

#
# Main program here
#

atoms = create_atoms_N2H4()
#println(atoms)

pw = PWGrid(20.0, atoms.LatVecs, dual=5.0)
#println(pw)

Nspecies = atoms.Nspecies    
pspfiles = [
    "/home/efefer/pseudo/PSLIB/N.pbe-n-rrkjus_psl.0.1.UPF",
    "/home/efefer/pseudo/PSLIB/H.pbe-rrkjus_psl.0.1.UPF"
]

pspots = Vector{PsPot_UPF}(undef,Nspecies)
for isp = 1:Nspecies
    pspots[isp] = PsPot_UPF( pspfiles[isp] )
    PWDFT._build_prj_interp_table!( pspots[isp], pw )
end


PsPotNL_UPF(atoms, pw, pspots)

