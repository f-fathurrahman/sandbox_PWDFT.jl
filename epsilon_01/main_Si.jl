using Printf
using PWDFT
import Serialization

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function init_Hamiltonian()    
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 20.0
    return Hamiltonian( atoms, pspfiles, ecutwfc,
        xcfunc="VWN", meshk=[6,6,6], shiftk=[1,1,1] ) # for FFT grid to have the same size as pwscf
end


function main_scf()
    Random.seed!(1234)
    Ham = init_Hamiltonian()
    println(Ham)
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.5 )
    Serialization.serialize("Ham.data", Ham)
end

function main_nscf()

    Random.seed!(1234)

    Ham = Serialization.deserialize("Ham.data")

    atoms = Ham.atoms
    ecutwfc = Ham.pw.ecutwfc
    Nspin = Ham.electrons.Nspin

    # more dense k-point grid, no symmetry
    kpoints = KPoints( atoms, [10,10,10], [1,1,1], SymmetryInfo().s, time_reversal=false )

    # New pw
    pw = PWGrid(ecutwfc, atoms.LatVecs, kpoints=kpoints)
    Ham.pw = pw

    Ngw = Ham.pw.gvecw.Ngw
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    k = Ham.pw.gvecw.kpoints.k

    # Manually construct Ham.electrons, add more empty states
    Ham.electrons = Electrons( atoms, Ham.pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                               Nstates_empty=16 )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    Nspin = Ham.electrons.Nspin

    Nkspin = Nkpt*Nspin

    # Also, new PsPotNL
    Ham.pspotNL = PsPotNL( atoms, pw, Ham.pspots )

    psiks = rand_BlochWavefunc( Ham )
    evals = zeros(Float64,Nstates,Nkspin)

    # Verbose
    k = Ham.pw.gvecw.kpoints.k
    for ispin = 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        #
        @printf("\nispin = %d, ik = %d, ikspin=%d, Ngw = %d\n", ispin, ik, i, Ngw[ik])
        @printf("kpts = [%f,%f,%f]\n", k[1,ik], k[2,ik], k[3,ik])
        evals[:,i] = diag_LOBPCG!( Ham, psiks[i], verbose_last=true )
    end

    Serialization.serialize("Ham_nscf.data", Ham)
    Serialization.serialize("psiks.data", psiks)
    Serialization.serialize("evals.data", evals)

end

main_scf()
main_nscf()
