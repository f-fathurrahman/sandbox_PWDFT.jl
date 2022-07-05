using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

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


function test_main()
    atoms, pw, pspots, electrons = init_from_pwinput()
    
    pspotNL = PsPotNL_UPF(atoms, pw, pspots)
    pspotNL_ref = PsPotNL(atoms, pw, pspots)

    # Check Vnl_KB construction
    ik = 1
    nkb = pspotNL.nkb
    Vnl_KB = zeros(ComplexF64, pw.gvecw.Ngw[ik], nkb)
    _init_Vnl_KB!( ik, atoms, pw, pspots, pspotNL, Vnl_KB )

    Ngwk = pw.gvecw.Ngw[ik]
    Nstates = electrons.Nstates

    psi = ones(ComplexF64, Ngwk, Nstates)

    println("\nUsing PsPotNL_UPF")
    println("-----------------")
    println("sum betaNL = ", sum(Vnl_KB))
    println("size betaNL = ", size(Vnl_KB))
    betaNL_psi = Vnl_KB' * psi
    println("sum betaNL_psi = ", sum(betaNL_psi))
    #display(abs.(betaNL_psi[1:5,1:5])); println();

    println("\nUsing PsPotNL")
    println("-------------")
    println("sum betaNL = ", sum(pspotNL_ref.betaNL[1]))
    println("size betaNL = ", size(pspotNL_ref.betaNL[1]))
    betaNL_psi = pspotNL_ref.betaNL[1]' * psi
    println("sum betaNL_psi = ", sum(betaNL_psi))
    #display(abs.(betaNL_psi[1:5,1:5])); println();

    println(pspotNL)

    println(pspotNL_ref)

end

test_main()
