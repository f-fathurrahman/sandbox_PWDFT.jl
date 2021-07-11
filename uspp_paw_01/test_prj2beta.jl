using Printf
using SpecialFunctions: sphericalbesselj
import LightXML

using PWDFT
const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("PsPot_UPF.jl")
include("integ_simpson.jl")

function test_prj2beta(
    atoms::Atoms, pspots
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m, lmax = 3 + 1

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m, lmax = 2 + 1

    prj2beta = Array{Int64}(undef,3,Natoms,4,7)
    prj2beta[:,:,:,:] .= -1   # set to invalid index

    NbetaNL = 0
    
    #for ia = 1:Natoms
    #    isp = atm2species[ia]
    #    psp = pspots[isp]
    #    for l in 0:psp.lmax
    #        for iprj in 1:psp.Nproj_l[l+1]
    #            #println("----")
    #            for m in -l:l
    #                NbetaNL = NbetaNL + 1
    #                prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
    #                @printf("ibeta=%3d ia=%3d l=%3d m=%3d\n", NbetaNL, ia, l, m)
    #            end
    #        end
    #    end
    #end

    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l in 0:psp.lmax, m in -l:l
            #println("------")
            for iprj in 1:psp.Nproj_l[l+1]
                NbetaNL = NbetaNL + 1
                prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
                @printf("ibeta=%3d ia=%3d l=%3d m=%3d\n", NbetaNL, ia, l, m)
            end
        end
    end

    println("NbetaNL = ", NbetaNL)

    return
end

function main_ZnO()
    atoms = Atoms( xyz_string_frac=
        """
        4

        Zn      0.3333333   0.6666667   0.0000000
        Zn      0.6666667   0.3333333   0.5000000
        O       0.3333333   0.6666667   0.3450000
        O       0.6666667   0.3333333   0.8450000
        """, in_bohr=true,
        LatVecs = gen_lattice_hexagonal( 3.2495*ANG2BOHR, 5.2069*ANG2BOHR ) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Zn-q2.gth"),
                joinpath(DIR_PSP, "O-q6.gth")]

    #pspfiles = ["/home/efefer/pseudo/HGH/Zn.pz-hgh.UPF",
    #            "/home/efefer/pseudo/HGH/O.pz-hgh.UPF"

    Nspecies = atoms.Nspecies
    @assert length(pspfiles) == Nspecies
    #pspots = Array{PsPot_UPF,1}(undef,Nspecies)
    pspots = Array{PsPot_GTH,1}(undef,Nspecies)
    for isp in 1:Nspecies
        pspots[isp] = PsPot_GTH(pspfiles[isp])
        #println(pspots[isp])
    end

    test_prj2beta(atoms, pspots)

end

main_ZnO()