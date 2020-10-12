using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("../get_default_psp.jl")

function test_O2()
    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", "O2.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    println(atoms)

    Nspecies = atoms.Nspecies
    pspots = Array{PsPot_GTH}(undef,Nspecies)
    for isp = 1:Nspecies
        pspots[isp] = PsPot_GTH(pspfiles[isp])
    end

    electrons = Electrons( atoms, pspots, (5,7) )
    println(electrons)

    electrons = Electrons( atoms, pspots, (7,5), Nstates_extra=1 )
    println(electrons)

end

test_O2()