using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

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

    # Check Vnl_KB construction
    ik = 1
    Vnl_KB = pspotNL.betaNL[ik]

    Ngwk = pw.gvecw.Ngw[ik]
    Nstates = electrons.Nstates
    
    psi = zeros(ComplexF64, Ngwk, Nstates)
    for i in 1:Nstates
        psi[i,i] = 1.0
    end

    betaNL_psi = Vnl_KB' * psi
    println("sum betaNL_psi = ", sum(betaNL_psi))

    display(abs.(betaNL_psi[1:5,1:5])); println();


    println("sum psi = ", sum(psi))

    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # Here the op_S operation
    nkb = pspotNL.NbetaNL
    nh = pspotNL.nh
    indv_ijkb0 = pspotNL.indv_ijkb0
    qq_at = pspotNL.qq_at

    ps = zeros(ComplexF64, nkb, Nstates)
    Spsi = copy(psi)

    for isp in 1:Nspecies
        #
        if pspots[isp].is_ultrasoft
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                #CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_DP,0.0_DP), &
                # qqc, nh(nt), becp%k(indv_ijkb0(na)+1,1), nkb, &
                # (0.0_DP,0.0_DP), ps(indv_ijkb0(na)+1,1), nkb )
                idx1 = (indv_ijkb0[ia]+1):(indv_ijkb0[ia]+nh[isp])
                ps[idx1,1:Nstates] = qq_at[1:nh[isp],1:nh[isp],ia] * betaNL_psi[idx1,:]

            end
        else
            if nh[isp] > 0
                for ia in 1:Natoms
                    if atm2species[ia] != isp
                        continue
                    end
                    idx1 = (indv_ijkb0[ia]+1):(indv_ijkb0[ia]+nh[isp])
                    ps[idx1,1:Nstates] .= 0.0 + im*0.0
                end
            end
        end
    end

    println("sum(ps) = ", sum(ps))

    Spsi = Spsi + Vnl_KB * ps

    println("sum Spsi = ", sum(Spsi))
    println("After: sum(abs(Spsi)) = ", sum(abs.(Spsi)) - sum(psi))

end

test_main()
