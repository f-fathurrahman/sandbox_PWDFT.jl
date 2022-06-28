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
function init_test_main()
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

    return atoms, pw, pspots
end

atoms, pw, pspots = init_test_main()
pspotNL = PsPotNL_UPF(atoms, pw, pspots)

println(pspotNL)

# Check Vnl_KB construction
ik = 1
nkb = pspotNL.nkb
Vnl_KB = zeros(ComplexF64, pw.gvecw.Ngw[ik], nkb)
_init_Vnl_KB!( ik, atoms, pw, pspots, pspotNL, Vnl_KB )

