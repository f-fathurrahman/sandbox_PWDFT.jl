using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("create_atoms_N2H4.jl")
include("usnldiag.jl")

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
    ispin = 1
    Vnl_KB = pspotNL.betaNL[ik]

    Ngwk = pw.gvecw.Ngw[ik]
    Nstates = electrons.Nstates
    psi = ones(ComplexF64, Ngwk, Nstates)
    betaNL_psi = Vnl_KB' * psi
    
    #if all(.!pspotNL.are_ultrasoft)
    #    println("No ultrasoft pspots")
    #    H_diag = zeros(Float64, Ngwk)
    #    usnldiag!( ik, ispin, atoms, pw, pspotNL, H_diag )
    #    println("sum H_diag = ", sum(H_diag))
    #else
    #    println("Using ultrasoft pspots")
        H_diag = zeros(Float64, Ngwk)
        S_diag = zeros(Float64, Ngwk)
        usnldiag!( ik, ispin, atoms, pw, pspotNL, H_diag, S_diag )
        println("sum H_diag = ", sum(H_diag))
        println("sum S_diag = ", sum(S_diag) - Ngwk)
    #end

end

test_main()
