using Printf
using Random
using LinearAlgebra
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("../get_default_psp.jl")
include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))
include("update_positions.jl")

Random.seed!(1234)

# From QE
const AMU_SI = 1.660538782e-27
const ELECTRONMASS_SI = 9.10938215e-31
const AMU_AU = AMU_SI / ELECTRONMASS_SI

const H_PLANCK_SI = 6.62607015e-34      # J s
const HARTREE_SI = 4.3597447222071e-18  # J
const AU_SEC = H_PLANCK_SI/(2*pi)/HARTREE_SI
const AU_PS  = AU_SEC*1.0e12

function init_Ham_H2O()
    ecutwfc = 15.0
    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", "H2O.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    # Set masses
    Ham.atoms.masses[:] = [16.0, 2.0]*AMU_AU

    return Ham
end


function init_Ham_CO2()
    ecutwfc = 15.0
    #filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", "CO2.xyz")
    filename = joinpath("../md_01/CO2.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    # Set masses
    Ham.atoms.masses[:] = [14.0, 16.0]*AMU_AU

    return Ham
end

function run_pwdft_jl!( Ham, psis; NiterMax=100, etot_conv_thr=1e-8 )
    KS_solve_Emin_PCG!( Ham, psis,
        skip_initial_diag=true, etot_conv_thr=etot_conv_thr,
        NiterMax=NiterMax
    )
    forces = calc_forces( Ham, psis )
    return sum(Ham.energies), forces
end


function display_orthonormality(s::String, psi::Matrix{ComplexF64})
    Nstates = size(psi,2)
    @views s11 = dot(psi[:,1], psi[:,1])
    @views s12 = dot(psi[:,1], psi[:,2])
    println()
    println(s)
    println("dot psi1,psi1 = ", s11)
    println("dot psi1,psi2 = ", s12)
    return
end


function main( init_func; fnametrj="TRAJ.xyz", fnameetot="ETOT.dat" )

    dt_fs = 1.0
    # Time step, in Ha atomic unit
    dt = dt_fs*10e-16/AU_SEC
    println("dt (au) = ", dt)
    println("dt (ps) = ", dt*AU_PS)

    Ham = init_func()

    Nspin = Ham.electrons.Nspin

    # XL-BOMD parameters
    κ = 1.82
    ω2 = κ/dt^2
    α = 0.018  # α*10^-3 = 18 -> α = 18000 ? or 18x10^3 = 0.018
    c0 = -6.0
    c1 = 14.0
    c2 = -8.0
    c3 = -3.0
    c4 = 4.0
    c5 = -1.0

    println(Ham.atoms.masses)

    Natoms = Ham.atoms.Natoms
    # Momenta
    p = zeros(Float64,3,Natoms)

    dr = zeros(Float64,3,Natoms)
    v = zeros(Float64,3,Natoms)
    vtilde = zeros(Float64,3,Natoms)
    m = Ham.atoms.masses
    atm2species = Ham.atoms.atm2species

    filetraj = open(fnametrj, "w")
    fileetot = open(fnameetot, "w")

    #FORCE_evAng = 2*Ry2eV/BOHR2ANG
    # XXX Xcrysden assumes the forces are in Ha/angstrom
    FORCE_evAng = 1.0/BOHR2ANG

    Nstates = Ham.electrons.Nstates
    O = zeros(ComplexF64,Nstates,Nstates)
    U = zeros(ComplexF64,Nstates,Nstates)
    C = zeros(ComplexF64,Nstates,Nstates)

    psis = rand_BlochWavefunc(Ham)
    psis_SC = deepcopy(psis) # initialize memory
    # Set to zeros, useful for check whether phi's already initialized with previous values.
    phi_m0 = zeros_BlochWavefunc(Ham) # after minimized
    phi_m1 = zeros_BlochWavefunc(Ham) # initialize memory
    phi_m2 = zeros_BlochWavefunc(Ham) # initialize memory
    phi_m3 = zeros_BlochWavefunc(Ham) # initialize memory
    phi_m4 = zeros_BlochWavefunc(Ham) # initialize memory
    phi_m5 = zeros_BlochWavefunc(Ham) # initialize memory

    energies, forces = run_pwdft_jl!(Ham, psis_SC)
    psis_SC0 = deepcopy(psis_SC) # initialize memory

    # Initial condition for XL-BOMD
    for i in 1:Nspin
        phi_m0[i][:,:] = psis_SC[i][:,:]
    end

    Etot = sum(energies)

    println("Forces = ")
    for ia in 1:Natoms
        isp = atm2species[ia]
        atsymb = Ham.atoms.SpeciesSymbols[isp]
        @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
            forces[1,ia], forces[2,ia], forces[3,ia])
    end

    Ekin_ions = 0.0  # assume initial velocities is zeroes
    Etot_conserved = Etot + Ekin_ions

    #
    # Start MD loop here
    #
    NiterMax = 5000
    iter_start_XL = 5

    # 0-th iteration, initial condition
    iter = 0
    @printf(filetraj, "%d  Etot_conserved = %18.10f\n\n", Natoms, Etot_conserved)
    @printf(fileetot, "%18.10f %18.10f %18.10f %18.10f\n",
            AU_PS*iter*dt, Etot_conserved, Etot, Ekin_ions)
    for ia in 1:Natoms
        isp = atm2species[ia]
        r = Ham.atoms.positions
        @printf(filetraj, "%3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                Ham.atoms.SpeciesSymbols[isp],
                r[1,ia]*BOHR2ANG, r[2,ia]*BOHR2ANG, r[3,ia]*BOHR2ANG,
                forces[1,ia]*FORCE_evAng, forces[2,ia]*FORCE_evAng, forces[3,ia]*FORCE_evAng)
    end
    flush(filetraj)
    flush(fileetot)

    Etot = sum(energies)

    for ia in 1:Natoms
        isp = atm2species[ia]
        for i in 1:3
            v[i,ia] = vtilde[i,ia] + 0.5*dt*forces[i,ia]/m[isp]
        end
    end
        
    @printf("\nMD Iter = %3d, Etot = %18.10f\n", iter, sum(energies))
    println("Forces = ")
    for ia in 1:Natoms
        isp = atm2species[ia]
        atsymb = Ham.atoms.SpeciesSymbols[isp]
        @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
            forces[1,ia], forces[2,ia], forces[3,ia])
    end


    for iter = 1:NiterMax

        # Update atomic positions
        for ia in 1:Natoms
            isp = atm2species[ia]
            for i in 1:3
                vtilde[i,ia] = v[i,ia] + 0.5*dt*forces[i,ia]/m[isp]
                dr[i,ia] = dt*vtilde[i,ia]
            end
        end
        update_positions!( Ham, dr )

        if iter > iter_start_XL

            println("XL-BOMD step, check norms:")
            display_orthonormality("phi_m5", phi_m5[1])
            display_orthonormality("phi_m4", phi_m4[1])
            display_orthonormality("phi_m3", phi_m3[1])
            display_orthonormality("phi_m2", phi_m2[1])
            display_orthonormality("phi_m1", phi_m1[1])
            display_orthonormality("phi_m0", phi_m0[1])

            # phi dynamics
            for i in 1:Nspin
                # Alignment
                @views O[:,:] = psis_SC[i]' * phi_m0[i]
                @views U[:,:] = inv(sqrt(O*O')) * O
                #
                @views psis[i][:,:] = 2*phi_m0[i] - phi_m1[i] +
                    κ * ( psis_SC[i]*U - phi_m0[i] ) +
                    α * (
                        c0*phi_m0[i] + c1*phi_m1[i] +
                        c2*phi_m2[i] + c3*phi_m3[i] +
                        c4*phi_m4[i] + c5*phi_m5[i]
                    )
            end

            display_orthonormality("psis", psis[1])

            #
            # Prepare initial guess for the electronic minimization
            #
            for i in 1:Nspin
                @views psis_SC[i][:,:] = psis[i][:,:]
                # psis_SC should be in orthonormal states before passed to KS solver
                @views C[:,:] = inv( sqrt(psis_SC[i]' * psis_SC[i]) )
                @views psis_SC[i][:,:] = psis_SC[i][:,:]*C
            end
            energies, forces[:] = run_pwdft_jl!(Ham, psis_SC)

        else
            # The usual BOMD
            # For input to electronic minimization
            # reuse psis_SC from previous iteration

            # save previous psis_SC
            #for i in 1:Nspin
            #    psis_SC0[i][:,:] = psis_SC[i][:,:]
            #end
            energies, forces[:] = run_pwdft_jl!(Ham, psis_SC, etot_conv_thr=1e-10)
            # Alignment, for initial state
            for i in 1:Nspin
                O[:,:] = psis_SC[i]' * psis_SC0[i]
                U[:,:] = inv(sqrt(O*O')) * O
                psis_SC[i][:,:] = psis_SC[i][:,:]*U
            end
        end
        
        for i in 1:Nspin
            phi_m5[i][:,:] = phi_m4[i][:,:]
            phi_m4[i][:,:] = phi_m3[i][:,:]
            phi_m3[i][:,:] = phi_m2[i][:,:]
            phi_m2[i][:,:] = phi_m1[i][:,:]
            phi_m1[i][:,:] = phi_m0[i][:,:]
        end

        if iter > iter_start_XL
            # Save old wavefunction, for electron dynamics in XL-BOMD
            for i in 1:Nspin
                phi_m0[i][:,:] = psis[i][:,:]
            end
        else
            # For initial conditions in XLBOMD
            for i in 1:Nspin
                phi_m0[i][:,:] = psis_SC[i][:,:]
            end
        end

        Etot = sum(energies)

        # Update velocities
        Ekin_ions = 0.0
        for ia in 1:Natoms
            isp = atm2species[ia]
            ptot2 = 0.0
            for i in 1:3
                v[i,ia] = vtilde[i,ia] + 0.5*dt*forces[i,ia]/m[isp]
                ptot2 = ptot2 + v[i,ia]^2
            end
            Ekin_ions = Ekin_ions + 0.5*m[isp]*ptot2
        end
        
        Etot_conserved = Etot + Ekin_ions

        @printf(filetraj, "%d  Etot_conserved = %18.10f\n\n", Natoms, Etot_conserved)
        @printf(fileetot, "%18.10f %18.10f %18.10f %18.10f\n",
                AU_PS*iter*dt, Etot_conserved, Etot, Ekin_ions)
        for ia in 1:Natoms
            isp = atm2species[ia]
            r = Ham.atoms.positions
            @printf(filetraj, "%3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                    Ham.atoms.SpeciesSymbols[isp],
                    r[1,ia]*BOHR2ANG, r[2,ia]*BOHR2ANG, r[3,ia]*BOHR2ANG,
                    forces[1,ia]*FORCE_evAng, forces[2,ia]*FORCE_evAng, forces[3,ia]*FORCE_evAng)
        end
        flush(filetraj)
        flush(fileetot)

        @printf("\nMD Iter = %3d, Etot = %18.10f\n", iter, sum(energies))
        println("Forces = ")
        for ia in 1:Natoms
            isp = atm2species[ia]
            atsymb = Ham.atoms.SpeciesSymbols[isp]
            @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
                forces[1,ia], forces[2,ia], forces[3,ia])
        end

    end

    close(filetraj)
    close(fileetot)

end

#main(init_Ham_H2O, fnametrj="TRAJ_H2O_v4.xyz", fnameetot="ETOT_H2O_v4.dat")
main(init_Ham_CO2,
    fnametrj="TRAJ_CO2_xlbomd.xyz",
    fnameetot="ETOT_CO2_xlbomd.dat"
)
#main(init_Ham_CO2,
#    fnametrj="TRAJ_CO2_step10_extrap2nd.xyz",
#    fnameetot="ETOT_CO2_step10_extrap2nd.dat"
#)
