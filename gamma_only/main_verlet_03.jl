using Printf
using Random
using LinearAlgebra
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))
include("INC_gamma_only.jl")
include("update_positions.jl")
include("find_extrap_alpha_beta.jl")

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
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )
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
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [14.0, 16.0]*AMU_AU

    return Ham
end

function run_pwdft_jl!( Ham, psis )
    KS_solve_Emin_PCG_dot!( Ham, psis, skip_initial_diag=true, etot_conv_thr=1e-8 )
    forces = calc_forces( Ham, psis )
    return sum(Ham.energies), forces
end

function main( init_func; fnametrj="TRAJ.xyz", fnameetot="ETOT.dat" )

    dt_fs = 1.0
    # Time step, in Ha atomic unit
    dt = dt_fs*10e-16/AU_SEC
    println("dt (au) = ", dt)
    println("dt (ps) = ", dt*AU_PS)

    Ham = init_func()

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

    tau_old = zeros(3,Natoms,3)
    tau = zeros(3,Natoms)

    psis = randn_BlochWavefuncGamma(Ham)
    energies, forces = run_pwdft_jl!(Ham, psis)
    
    psis_m0 = deepcopy(psis) # after minimized
    psis_m1 = deepcopy(psis) # initialize memory
    psis_m2 = deepcopy(psis) # initialize memory

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

    Nspin = Ham.electrons.Nspin

    Nstates = Ham.electrons.Nstates
    O = zeros(ComplexF64,Nstates,Nstates)
    U = zeros(ComplexF64,Nstates,Nstates)

    #
    # Start MD loop here
    #
    NiterMax = 20
    for iter = 1:NiterMax

        @printf(filetraj, "%d  Etot_conserved = %18.10f\n\n", Natoms, Etot_conserved)
        @printf(fileetot, "%18.10f %18.10f %18.10f %18.10f\n",
            AU_PS*(iter-1)*dt, Etot_conserved, Etot, Ekin_ions)
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

        tau_old[:,:,3] = tau_old[:,:,2]
        tau_old[:,:,2] = tau_old[:,:,1]
        tau_old[:,:,1] = Ham.atoms.positions[:,:]

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
        Etot_conserved = Etot + Ekin_ions
        update_positions!( Ham, Ham.atoms.positions + dr )

        tau[:,:] = Ham.atoms.positions[:,:]

        # Extrapolate wavefunction (1st order)
        if iter > 3
            α, β = find_extrap_alpha_beta(tau, tau_old)
            for i in 1:Nspin
                psis.data[i][:,:] = psis_m0.data[i][:,:] +
                                    α*(psis_m0.data[i][:,:] - psis_m1.data[i][:,:]) +
                                    β*(psis_m1.data[i][:,:] - psis_m2.data[i][:,:])
                ortho_sqrt_gamma!( psis.data[i] )
            end
            #println("Check ortho: ", dot(psis,psis))
        end

        for i in 1:Nspin
            psis_m2.data[i][:,:] = psis_m1.data[i][:,:]
            psis_m1.data[i][:,:] = psis_m0.data[i][:,:]
        end
        
        energies, forces[:] = run_pwdft_jl!(Ham, psis)
        
        # Alignment
        for i in 1:Nspin
            O[:,:] = overlap_gamma(psis.data[i], psis_m0.data[i])
            U[:,:] = inv(sqrt(O' * O)) * O'
            psis.data[i] = psis.data[i]*U
        end

        for i in 1:Nspin
            psis_m0.data[i][:,:] = psis.data[i][:,:]
        end

        Etot = sum(energies)

        for ia in 1:Natoms
            isp = atm2species[ia]
            for i in 1:3
                v[i,ia] = vtilde[i,ia] + 0.5*dt*forces[i,ia]/m[isp]
            end
        end
        
        @printf("\nIter = %3d, Etot = %18.10f\n", iter, sum(energies))
        println("Forces = ")
        for ia in 1:Natoms
            isp = atm2species[ia]
            atsymb = Ham.atoms.SpeciesSymbols[isp]
            @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
                forces[1,ia], forces[2,ia], forces[3,ia])
        end

    end

    @printf(filetraj, "%d  Etot_conserved = %18.10f\n\n", Natoms, Etot_conserved)
    @printf(fileetot, "%18.10f %18.10f %18.10f %18.10f\n",
        AU_PS*NiterMax*dt, Etot_conserved, Etot, Ekin_ions)
    for ia in 1:Natoms
        isp = atm2species[ia]
        r = Ham.atoms.positions
        @printf(filetraj, "%3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                Ham.atoms.SpeciesSymbols[isp],
                r[1,ia]*BOHR2ANG, r[2,ia]*BOHR2ANG, r[3,ia]*BOHR2ANG,
                forces[1,ia]*FORCE_evAng, forces[2,ia]*FORCE_evAng, forces[3,ia]*FORCE_evAng)
    end

    close(filetraj)
    close(fileetot)

end

#main(init_Ham_H2O, fnametrj="TRAJ_H2O_v4.xyz", fnameetot="ETOT_H2O_v4.dat")
main(init_Ham_CO2,
    fnametrj="TRAJ_CO2_extrap2nd.xyz",
    fnameetot="ETOT_CO2_extrap2nd.dat"
)
#main(init_Ham_CO2,
#    fnametrj="TRAJ_CO2_step10_extrap2nd.xyz",
#    fnameetot="ETOT_CO2_step10_extrap2nd.dat"
#)