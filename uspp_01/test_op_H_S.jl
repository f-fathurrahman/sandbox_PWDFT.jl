using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("update_from_rhoe.jl")
include("newd.jl")

function op_S!(
    Ham::Hamiltonian{Txc,PsPot_UPF},
    psi::AbstractArray{ComplexF64},
    Spsi::AbstractArray{ComplexF64}
) where Txc <: AbstractXCCalculator

    # Check Vnl_KB construction
    ik = Ham.ik
    ispin = Ham.ispin # XXX: used?

    pspotNL = Ham.pspotNL
    pspots = Ham.pspots
    atoms = Ham.atoms

    Ngwk = size(psi,1) # check againsts Ngw[k]
    Nstates = size(psi,2)

    Vnl_KB = pspotNL.betaNL[ik]

    betaNL_psi = Vnl_KB' * psi  # XXX precalculate this
    
    #println("sum betaNL_psi = ", sum(betaNL_psi))
    #display(abs.(betaNL_psi[1:5,1:5])); println();

    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nkb = pspotNL.NbetaNL
    nh = pspotNL.nh
    indv_ijkb0 = pspotNL.indv_ijkb0
    qq_at = pspotNL.qq_at

    ps = zeros(ComplexF64, nkb, Nstates)
    
    @views Spsi[:,:] = psi[:,:]

    for isp in 1:Nspecies
        #
        if pspots[isp].is_ultrasoft
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
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

    @views Spsi[:,:] += Vnl_KB * ps

    return

end


#atoms, pw, pspots = init_test_main()
function test_main()

    Ham = init_Ham_from_pwinput()


    Rhoe, RhoeG = atomic_rho_g(Ham)

    update_from_rhoe!( Ham, Rhoe, RhoeG )

    calc_newDeeq!( Ham )

    ik = 1
    Ham.ik = ik # set the current k index of H
    Ngwk = Ham.pw.gvecw.Ngw[ik]
    Nstates = Ham.electrons.Nstates

    psi = zeros(ComplexF64, Ngwk, Nstates)
    for i in 1:Nstates
        psi[i,i] = 1.0
    end

    Spsi = zeros(ComplexF64, Ngwk, Nstates)
    Hpsi = zeros(ComplexF64, Ngwk, Nstates)

    op_H!(Ham, psi, Hpsi)
    op_S!(Ham, psi, Spsi)

    println("sum Hpsi = ", sum(Hpsi))
    println("sum(abs(Hpsi)) = ", sum(abs.(Hpsi)))

    println("sum Spsi = ", sum(Spsi))
    println("sum(abs(Spsi)) min ref = ", sum(abs.(Spsi)) - sum(psi))

end

test_main()
