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

Random.seed!(1234)

# From QE
const AMU_SI = 1.660538782e-27  # kg
const ELECTRONMASS_SI = 9.10938215e-31 # kg
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

function init_Ham_Si8()
    ecutwfc = 15.0
    atoms = Atoms(xyz_string="""
    8

    Si       0.0000000000       0.0000000000       0.0000000000
    Si       7.6969017521       7.6969017521       2.5656339174
    Si       5.1312678348       0.0000000000       5.1312678348
    Si       7.6969017521       2.5656339174       7.6969017521
    Si       0.0000000000       5.1312678348       5.1312678348
    Si       2.5656339174       2.5656339174       2.5656339174
    Si       2.5656339174       7.6969017521       7.6969017521
    Si       5.1312678348       5.1312678348       0.0000000000
    """, in_bohr=true, LatVecs=gen_lattice_sc(10.2625356695))
    write_xsf("TEMP_Si8.xsf", atoms)
    pspfiles = get_default_psp(atoms)
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [28.085]*AMU_AU
    return Ham
end


function init_Ham_Si64()
    ecutwfc = 15.0
    atoms = Atoms(xyz_string_frac="""
    64

    Si   0.000000000000000   0.000000000000000   0.000000000000000 
    Si   0.375000000000000   0.375000000000000   0.125000000000000 
    Si   0.250000000000000   0.000000000000000   0.250000000000000 
    Si   0.375000000000000   0.125000000000000   0.375000000000000 
    Si   0.000000000000000   0.250000000000000   0.250000000000000 
    Si   0.125000000000000   0.125000000000000   0.125000000000000 
    Si   0.125000000000000   0.375000000000000   0.375000000000000 
    Si   0.250000000000000   0.250000000000000   0.000000000000000 
    Si   0.500000000000000   0.500000000000000   0.000000000000000 
    Si   0.500000000000000   0.000000000000000   0.000000000000000 
    Si   0.000000000000000   0.500000000000000   0.500000000000000 
    Si   0.000000000000000   0.000000000000000   0.500000000000000 
    Si   0.500000000000000   0.500000000000000   0.500000000000000 
    Si   0.500000000000000   0.000000000000000   0.500000000000000 
    Si   0.000000000000000   0.500000000000000   0.000000000000000 
    Si   0.875000000000000   0.875000000000000   0.125000000000000 
    Si   0.875000000000000   0.375000000000000   0.125000000000000 
    Si   0.375000000000000   0.875000000000000   0.625000000000000 
    Si   0.375000000000000   0.375000000000000   0.625000000000000 
    Si   0.875000000000000   0.875000000000000   0.625000000000000 
    Si   0.875000000000000   0.375000000000000   0.625000000000000 
    Si   0.375000000000000   0.875000000000000   0.125000000000000 
    Si   0.750000000000000   0.500000000000000   0.250000000000000 
    Si   0.750000000000000   0.000000000000000   0.250000000000000 
    Si   0.250000000000000   0.500000000000000   0.750000000000000 
    Si   0.250000000000000   0.000000000000000   0.750000000000000 
    Si   0.750000000000000   0.500000000000000   0.750000000000000 
    Si   0.750000000000000   0.000000000000000   0.750000000000000 
    Si   0.250000000000000   0.500000000000000   0.250000000000000 
    Si   0.875000000000000   0.625000000000000   0.375000000000000 
    Si   0.875000000000000   0.125000000000000   0.375000000000000 
    Si   0.375000000000000   0.625000000000000   0.875000000000000 
    Si   0.375000000000000   0.125000000000000   0.875000000000000 
    Si   0.875000000000000   0.625000000000000   0.875000000000000 
    Si   0.875000000000000   0.125000000000000   0.875000000000000 
    Si   0.375000000000000   0.625000000000000   0.375000000000000 
    Si   0.500000000000000   0.750000000000000   0.250000000000000 
    Si   0.500000000000000   0.250000000000000   0.250000000000000 
    Si   0.000000000000000   0.750000000000000   0.750000000000000 
    Si   0.000000000000000   0.250000000000000   0.750000000000000 
    Si   0.500000000000000   0.750000000000000   0.750000000000000 
    Si   0.500000000000000   0.250000000000000   0.750000000000000 
    Si   0.000000000000000   0.750000000000000   0.250000000000000 
    Si   0.625000000000000   0.625000000000000   0.125000000000000 
    Si   0.625000000000000   0.125000000000000   0.125000000000000 
    Si   0.125000000000000   0.625000000000000   0.625000000000000 
    Si   0.125000000000000   0.125000000000000   0.625000000000000 
    Si   0.625000000000000   0.625000000000000   0.625000000000000 
    Si   0.625000000000000   0.125000000000000   0.625000000000000 
    Si   0.125000000000000   0.625000000000000   0.125000000000000 
    Si   0.625000000000000   0.875000000000000   0.375000000000000 
    Si   0.625000000000000   0.375000000000000   0.375000000000000 
    Si   0.125000000000000   0.875000000000000   0.875000000000000 
    Si   0.125000000000000   0.375000000000000   0.875000000000000 
    Si   0.625000000000000   0.875000000000000   0.875000000000000 
    Si   0.625000000000000   0.375000000000000   0.875000000000000 
    Si   0.125000000000000   0.875000000000000   0.375000000000000 
    Si   0.750000000000000   0.750000000000000   0.000000000000000 
    Si   0.750000000000000   0.250000000000000   0.000000000000000 
    Si   0.250000000000000   0.750000000000000   0.500000000000000 
    Si   0.250000000000000   0.250000000000000   0.500000000000000 
    Si   0.750000000000000   0.750000000000000   0.500000000000000 
    Si   0.750000000000000   0.250000000000000   0.500000000000000 
    Si   0.250000000000000   0.750000000000000   0.000000000000000""",
    in_bohr=true, LatVecs=gen_lattice_sc(2*10.2625356695))
    write_xsf("TEMP_Si64.xsf", atoms)
    pspfiles = get_default_psp(atoms)
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [28.085]*AMU_AU
    return Ham
end


function init_Ham_C8()
    ecutwfc = 15.0
    atoms = Atoms(xyz_string_frac="""
    8

    C   0.000000000000000   0.000000000000000   0.000000000000000 
    C   0.750000000000000   0.750000000000000   0.250000000000000 
    C   0.500000000000000   0.000000000000000   0.500000000000000 
    C   0.750000000000000   0.250000000000000   0.750000000000000 
    C   0.000000000000000   0.500000000000000   0.500000000000000 
    C   0.250000000000000   0.250000000000000   0.250000000000000 
    C   0.250000000000000   0.750000000000000   0.750000000000000 
    C   0.500000000000000   0.500000000000000   0.000000000000000 
    """, in_bohr=true, LatVecs=gen_lattice_sc(3.57*ANG2BOHR))
    write_xsf("TEMP_C8.xsf", atoms)
    pspfiles = get_default_psp(atoms)
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )
    # Set masses
    Ham.atoms.masses[:] = [28.085]*AMU_AU
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

    # t = 0 for MD
    psis = randn_BlochWavefuncGamma(Ham)
    energies, forces = run_pwdft_jl!(Ham, psis)
    Etot = sum(energies)

    println("Forces = ")
    for ia in 1:Natoms
        isp = atm2species[ia]
        atsymb = Ham.atoms.SpeciesSymbols[isp]
        @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
            forces[1,ia], forces[2,ia], forces[3,ia])
    end

    Ekin_ions = 0.0  # assume initial velocities are zeroes
    Etot_conserved = Etot + Ekin_ions

    #
    # Start MD loop here
    #
    NiterMax = 5000
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

        energies, forces[:] = run_pwdft_jl!(Ham, psis)
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
#main(init_Ham_CO2,
#    fnametrj="TRAJ_CO2_noextrap.xyz",
#    fnameetot="ETOT_CO2_noextrap.dat"
#)

#main(init_Ham_C8,
#    fnametrj="TRAJ_C8_noextrap.xyz",
#    fnameetot="ETOT_C8_noextrap.dat"
#)

main(init_Ham_Si64,
    fnametrj="TRAJ_Si64_noextrap.xyz",
    fnameetot="ETOT_Si64_noextrap.dat"
)

#main(init_Ham_CO2,
#    fnametrj="TRAJ_CO2_step10.xyz",
#    fnameetot="ETOT_CO2_step10.dat"
#)