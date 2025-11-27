using Infiltrator
using Printf
using LinearAlgebra
using Random
using PWDFT
import Serialization

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

include(joinpath(DIR_PWDFT, "examples", "common", "dump_bandstructure.jl"))
include("my_gen_kpath.jl")

function main_bands()

    Ham = Serialization.deserialize("Ham.data")

    atoms = Ham.atoms
    ecutwfc = Ham.pw.ecutwfc

    # Band structure calculation
    kpoints, kpt_spec, kpt_spec_labels = my_gen_kpath(atoms, "G-S-N-S1-Z"; Δk=0.02)
    println("Nkpt = ", kpoints.Nkpt)

    # New pw
    pw = PWGrid(ecutwfc, atoms.LatVecs, kpoints=kpoints)
    Ham.pw = pw

    Ngw = Ham.pw.gvecw.Ngw
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    k = Ham.pw.gvecw.kpoints.k

    # Manually construct Ham.electrons
    Ham.electrons = Electrons( atoms, Ham.pspots, Nspin=1, Nkpt=kpoints.Nkpt,
                               Nstates_empty=5 )

    ebands = Ham.electrons.ebands
    Nspin = Ham.electrons.Nspin

    # Also, new PsPotNL
    Ham.pspotNL = PsPotNL( atoms, pw, Ham.pspots )

    psiks = rand_BlochWavefunc( Ham )

    k = Ham.pw.gvecw.kpoints.k
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        #
        @printf("\nispin = %d, ik = %d, ikspin=%d, Ngw = %d\n", ispin, ik, ikspin, Ngw[ik])
        @printf("kpts = [%f,%f,%f]\n", k[1,ik], k[2,ik], k[3,ik])
        ebands[:,ikspin], psiks[ikspin] =
        diag_LOBPCG( Ham, psiks[ikspin], verbose_last=true )
    end

    #@infiltrate
    dump_bandstructure( ebands, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_bands.dat" )

end

#main_bands()
