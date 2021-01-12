using Printf
using LinearAlgebra
using Random
using FFTW
using LightXML

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures/DATA_G2_mols")

include("PsPot_UPF.jl")
include("PsPotNL_UPF.jl")

function test_main()
    Random.seed!(1234)
    filename = joinpath(DIR_STRUCTURES, "N2H4.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspots = Array{PsPot_UPF,1}(undef, atoms.Nspecies)
    #pspots[1] = PsPot_UPF("GBRV_LDA/n_lda_v1.2.uspp.F.UPF2")
    #pspots[2] = PsPot_UPF("GBRV_LDA/h_lda_v1.4.uspp.F.UPF2")

    pspots[1] = PsPot_UPF("/home/efefer/pseudo/HGH/N.pz-hgh.UPF")
    pspots[2] = PsPot_UPF("/home/efefer/pseudo/HGH/H.pz-hgh.UPF")

    print(pspots[1])
    print(pspots[2])

    PsPotNL_UPF(atoms, pspots)
end

test_main()

