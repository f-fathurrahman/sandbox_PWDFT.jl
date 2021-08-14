using Printf
using PWDFT
import Serialization

include("calc_dipole_matrix.jl")

function load_Hamiltonian()
    return Serialization.deserialize("Ham.data")
end

function create_kpt_grid(atoms::Atoms; kgrid=[10,10,10])
    # more dense k-point grid, no symmetry
    kpoints = KPoints( atoms, kgrid, [1,1,1],
        SymmetryInfo().s, time_reversal=false
    )
    return kpoints
end


function update_Ham_k!(Ham, ik, kpoints_dense, Nstates_empty)

    atoms = Ham.atoms
    ecutwfc = Ham.pw.ecutwfc
    gvec = Ham.pw.gvec
    Nspin = Ham.electrons.Nspin

    # Prepare one-kpoint from kpoints_dense
    k = zeros(3,1)
    k[1:3,1] = kpoints_dense.k[:,ik]
    RecVecs = Ham.pw.RecVecs
    kpoints = KPoints(1, (0,0,0), k, [1.0], RecVecs)

    # New gvecw
    # FIXME: does not work currently, as PWGrid is immutable
    #Ham.pw.gvecw = PWDFT.init_gvecw( ecutwfc, gvec, kpoints )

    # New pw
    Ham.pw = PWGrid(ecutwfc, atoms.LatVecs, kpoints=kpoints)

    Ngw = Ham.pw.gvecw.Ngw
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    k = Ham.pw.gvecw.kpoints.k

    # Manually construct Ham.electrons, add more empty states
    Ham.electrons = Electrons( atoms, Ham.pspots, Nspin=Nspin, Nkpt=1,
                               Nstates_empty=16 )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    Nspin = Ham.electrons.Nspin

    Nkspin = Nkpt*Nspin

    # Also, new PsPotNL
    Ham.pspotNL = PsPotNL( atoms, Ham.pw, Ham.pspots )

    return
end

function main_nscf_dipole()

    Ham = load_Hamiltonian()

    kpoints_dense = create_kpt_grid( Ham.atoms, kgrid=[15,15,15] )

    Serialization.serialize("atoms.data", Ham.atoms) # for CellVolume
    Serialization.serialize("kpoints_nscf.data", kpoints_dense)

    Nstates_empty = 16 # XXX
    Nkpt_dense = kpoints_dense.Nkpt

    # We calculate new total states here
    Nstates = Ham.electrons.Nstates_occ + Nstates_empty

    evals = zeros(Float64,Nstates,Nkpt_dense)

    Ngwx = round(Int64,Ham.pw.gvec.Ng/6) # FIXME: properly calculate Ngwx
    psi_buffer = rand(ComplexF64,Ngwx,Nstates)

    M_aux = zeros(ComplexF64,3,Nstates,Nstates)
    Mk = Array{Array{Float64,3},1}(undef,Nkpt_dense)
    for ik in 1:Nkpt_dense
        Mk[ik] = zeros(Float64,3,Nstates,Nstates)
    end

    for ik in 1:Nkpt_dense
        println()
        #
        print("Update Ham k: ")
        @time update_Ham_k!(Ham, ik, kpoints_dense, Nstates_empty)
        #
        Npw = Ham.pw.gvecw.Ngw[1]
        @views psi = psi_buffer[1:Npw,:]
        #
        Ham.ik = 1
        Ham.ispin = 1
        print("Diagonalization: ")
        @time evals[:,ik] = diag_LOBPCG!( Ham, psi, verbose_last=false )
        #
        @printf("Done: ik = %8d, Npw = %8d Ngwx = %8d\n", ik, Npw, Ngwx)
        #
        print("Calculating dipole matrix: ")
        M = Mk[ik]
        @time begin
            calc_dipole_matrix!(
                Ham.pw, psi,
                1, Ham.electrons.Focc, ik, M_aux; metal_like=false
            )
        end
        for i in 1:length(M)
            M[i] = real( M_aux[i] * conj(M_aux[i]) )
        end
    end

    Serialization.serialize("Mk.data", Mk)
    Serialization.serialize("evals_nscf.data", evals)
    Serialization.serialize("electrons_nscf.data", Ham.electrons)

end

main_nscf_dipole()
