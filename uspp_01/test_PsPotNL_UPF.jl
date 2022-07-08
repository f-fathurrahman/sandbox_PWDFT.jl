using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("create_atoms_N2H4.jl")
include("../ylm_real/Ylm_real_qe.jl")
include("calc_clebsch_gordan.jl")
include("calc_qradG.jl")
include("qvan2.jl")
include("PsPotNL_UPF.jl")

include("../pwscf_02/PWSCFInput.jl")

function init_from_pwinput()
    println("ARGS = ", ARGS)
    @assert length(ARGS) == 1
    pwinput = PWSCFInput(ARGS[1])

    atoms = pwinput.atoms
    println(atoms)

    dual = pwinput.ecutrho/pwinput.ecutwfc
    pw = PWGrid(pwinput.ecutwfc, atoms.LatVecs, dual=dual)
    println(pw)

    Nspecies = atoms.Nspecies    
    pspfiles = pwinput.pspfiles
    pspots = Vector{PsPot_UPF}(undef,Nspecies)
    for isp in 1:Nspecies
        pspots[isp] = PsPot_UPF( pspfiles[isp] )
        PWDFT._build_prj_interp_table!( pspots[isp], pw )
    end

    electrons = Electrons(
        atoms, pspots,
        Nspin=1, Nkpt=pw.gvecw.kpoints.Nkpt, Nstates_empty=0
    )

    return atoms, pw, pspots, electrons
end



#atoms, pw, pspots = init_test_main()
function test_main()
    atoms, pw, pspots, electrons = init_from_pwinput()
    pspotNL = PsPotNL_UPF(atoms, pw, pspots)

    println(pspotNL)

    # Check Vnl_KB construction
    ik = 1
    Vnl_KB = pspotNL.betaNL[ik]

    Ngwk = pw.gvecw.Ngw[ik]
    Nstates = electrons.Nstates
    psi = ones(ComplexF64, Ngwk, Nstates)
    betaNL_psi = Vnl_KB' * psi
    println("sum betaNL_psi = ", sum(betaNL_psi))

    #display(abs.(betaNL_psi[1:5,1:5])); println();

end

test_main()
