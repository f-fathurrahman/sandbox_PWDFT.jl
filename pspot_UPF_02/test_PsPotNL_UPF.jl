using Printf
using LinearAlgebra
using Random
using FFTW
using LightXML

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures/DATA_G2_mols")

const UPF_HGH_DIR = "/home/efefer/pseudo/HGH"

include("PsPot_UPF.jl")
include("PsPotNL_UPF.jl")

function create_pspots_UPF_HGH(atoms::Atoms; xcname="pz")
    Nspecies = atoms.Nspecies
    SpeciesSymbols = atoms.SpeciesSymbols
    pspots = Array{PsPot_UPF}(undef,Nspecies)
    for isp in 1:Nspecies
        filename = SpeciesSymbols[isp]*"."*xcname*"-hgh.UPF"
        pspots[isp] = PsPot_UPF(joinpath(UPF_HGH_DIR, filename))
    end
    return pspots
end

function test_main()
    #filename = joinpath(DIR_STRUCTURES, "N2H4.xyz")
    #filename = joinpath(DIR_STRUCTURES, "N2.xyz")
    filename = joinpath(DIR_STRUCTURES, "Cl2.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspots = Array{PsPot_UPF,1}(undef, atoms.Nspecies)
    #pspots[1] = PsPot_UPF("GBRV_LDA/n_lda_v1.2.uspp.F.UPF2")
    #pspots[2] = PsPot_UPF("GBRV_LDA/h_lda_v1.4.uspp.F.UPF2")

    pspots = create_pspots_UPF_HGH(atoms)

    for pspot in pspots
        print(pspot)
    end
    PsPotNL_UPF(atoms, pspots)
end

test_main()

