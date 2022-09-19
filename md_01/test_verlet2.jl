using Printf
using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

include("../NLopt_v3/calc_energies_grad.jl")
include("../NLopt_v3/KS_solve_Emin_PCG_new.jl")
include("../NLopt_v3/KS_solve_Emin_PCG_vec.jl")
include("update_positions.jl")

Random.seed!(1234)

# From QE
const AMU_SI = 1.660538782e-27
const ELECTRONMASS_SI = 9.10938215e-31
const AMU_AU = AMU_SI / ELECTRONMASS_SI

#FORCE_evAng = 2*Ry2eV/BOHR2ANG
# XXX Xcrysden assumes the forces are in Ha/angstrom
const FORCE_evAng = 1.0/BOHR2ANG


function init_Ham_H2O()
    # Atoms
    atoms = Atoms( ext_xyz_file="H2O.xyz" )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "O-q6.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [16.0, 2.0]*AMU_AU

    return Ham
end


function init_Ham_CO2()
    # Atoms
    atoms = Atoms( ext_xyz_file="CO2.xyz" )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "C-q4.gth"),
                joinpath(DIR_PSP, "O-q6.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    # Set masses
    Ham.atoms.masses[:] = [14.0, 16.0]*AMU_AU

    return Ham
end

function run_pwdft_jl( Ham )
    psiks = rand_BlochWavefunc(Ham)
    #KS_solve_Emin_PCG!( Ham, psiks )
    KS_solve_Emin_PCG_vec!( Ham, psiks )
    forces = calc_forces( Ham, psiks )
    return sum(Ham.energies), forces
end

function run_pwdft_jl!( Ham, psiks )
    KS_solve_Emin_PCG_vec!( Ham, psiks, skip_initial_diag=true, etot_conv_thr=1e-8 )
    #KS_solve_SCF_potmix!( Ham, psiks, etot_conv_thr=1e-8 ) # not good, slow convergence
    forces = calc_forces( Ham, psiks )
    return sum(Ham.energies), forces
end

function main( init_func; fnametrj="TRAJ.xyz", fnameetot="ETOT.dat" )
    
    Ham = init_func()

    println(Ham.atoms.masses)

    Natoms = Ham.atoms.Natoms

    psiks = rand_BlochWavefunc(Ham)
    energies, forces = run_pwdft_jl!(Ham, psiks)
    Etot = sum(energies)

    println("Initial r  =")
    display(Ham.atoms.positions'); println()
    println("Initial forces = ")
    display(forces'); println()

    # Momenta
    p = zeros(Float64,3,Natoms)
    Ekin_ions = 0.0  # assume initial velocities is zeroes

    dr = zeros(Float64,3,Natoms)
    v = zeros(Float64,3,Natoms)
    vtilde = zeros(Float64,3,Natoms)
    m = Ham.atoms.masses
    atm2species = Ham.atoms.atm2species

    # Time step
    dt = 10.0

    filetraj = open(fnametrj, "w")
    fileetot = open(fnameetot, "w")

    #
    # Start MD loop here
    #
    NiterMax = 10
    for iter in 1:NiterMax

        write_Etot_traj( filetraj, fileetot, (iter-1)*dt,
            Ham.atoms, Etot, Ekin_ions, forces,
        )

        Ekin_ions = 0.0
        for ia in 1:Natoms
            isp = atm2species[ia]
            ptot2 = 0.0
            for i in 1:3
                vtilde[i,ia] = v[i,ia] + 0.5*dt*forces[i,ia]/m[isp]
                dr[i,ia] = dt*vtilde[i,ia]
                ptot2 = ptot2 + v[i,ia]^2
            end
            Ekin_ions = Ekin_ions + 0.5*m[isp]*ptot2
        end

        # Update the positions and also the Hamiltonian
        update_positions!( Ham, dr )
        # Calculate the forces
        energies, forces[:] = run_pwdft_jl!(Ham, psiks)
        Etot = sum(energies)

        for ia in 1:Natoms
            isp = atm2species[ia]
            for i in 1:3
                v[i,ia] = vtilde[i,ia] + 0.5*dt*forces[i,ia]/m[isp]
            end
        end
        
        @printf("\nIter = %3d, Etot = %18.10f\n", iter, sum(energies))

    end

    write_Etot_traj( filetraj, fileetot, NiterMax*dt,
        Ham.atoms, Etot, Ekin_ions, forces )

    # Dont't forget to close the streams
    close(filetraj)
    close(fileetot)

    return
end



function write_Etot_traj(
    filetraj, fileetot, t::Float64,
    atoms, Etot, Ekin_ions, forces
)
    R = atoms.positions * BOHR2ANG
    atm2species = atoms.atm2species
    F = forces * FORCE_evAng
    spsymb = atoms.SpeciesSymbols
    Natoms = atoms.Natoms
    Etot_conserved = Etot + Ekin_ions
    #
    @printf(filetraj, "%d  Etot_conserved = %18.10f\n\n", Natoms, Etot_conserved)
    @printf(fileetot, "%18.10f %18.10f %18.10f %18.10f\n", t, Etot_conserved, Etot, Ekin_ions)
    #
    for ia in 1:Natoms
        isp = atm2species[ia]
        @printf(filetraj,
                "%3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                spsymb[isp],
                R[1,ia], R[2,ia], R[3,ia],
                F[1,ia], F[2,ia], F[3,ia])
    end
    # Immediately flush
    flush(filetraj)
    flush(fileetot)
    return
end

#main(init_Ham_H2O, fnametrj="TRAJ_H2O_v4.xyz", fnameetot="ETOT_H2O_v4.dat")
main(init_Ham_CO2, fnametrj="TRAJ_CO2_v3.xyz", fnameetot="ETOT_CO2_v3.dat")