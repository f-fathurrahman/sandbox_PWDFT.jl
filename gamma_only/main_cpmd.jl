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


function init_Ham_Si8()
    ecutwfc = 15.0
    atoms = Atoms(xyz_string="""
    8

    Si       0.1000000000       0.0000000000       0.0000000000
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


# Initial minimization
function minimize_electrons!( Ham, psis; NiterMax=200, etot_conv_thr=1e-8 )
    KS_solve_Emin_PCG_dot!( Ham, psis,
        skip_initial_diag=true, etot_conv_thr=etot_conv_thr,
        NiterMax=NiterMax
    )
    forces = calc_forces( Ham, psis )
    return sum(Ham.energies), forces
end

function calc_elec_energies!(Ham, psis, μ, dpsis)
    Rhoe = calc_rhoe(Ham, psis)
    update!(Ham, psis)
    energies = calc_energies(Ham, psis)
    Etot = sum(energies)
    Ekin_elec = 0.5*μ*dot(dpsis,dpsis)
    return Etot, Ekin_elec
end 


function write_traj_files(
    filetraj, fileetot,
    atoms, forces,
    Etot_conserved, Etot, Ekin_ions, Ekin_elec,
    dt, iter
)
    # Forces are in Ha/bohr
    # Xcrysden expect the forces are in Ha/angstrom
    F = forces*(1.0/BOHR2ANG)
    #
    @printf(filetraj, "%d  Etot_conserved = %18.10f\n\n", atoms.Natoms, Etot_conserved)
    @printf(fileetot, "%18.10f %18.10f %18.10f %18.10f %18.10f\n",
            AU_PS*iter*dt, Etot_conserved, Etot, Ekin_ions, Ekin_elec)
    for ia in 1:atoms.Natoms
        isp = atoms.atm2species[ia]
        r = atoms.positions*BOHR2ANG # convert from bohr to angstrom
        @printf(filetraj, "%3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                atoms.SpeciesSymbols[isp],
                r[1,ia], r[2,ia], r[3,ia],
                F[1,ia], F[2,ia], F[3,ia]
        ) # printf
    end
    flush(filetraj)
    flush(fileetot)

    return
end


function info_md(atoms, energies, forces, iter)
    @printf("\nMD Iter = %3d, Etot = %18.10f\n", iter, sum(energies))
    println("Forces = ")
    for ia in 1:atoms.Natoms
        isp = atoms.atm2species[ia]
        atsymb = atoms.SpeciesSymbols[isp]
        @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
            forces[1,ia], forces[2,ia], forces[3,ia])
    end
end


function calc_X_matrix(Ctilde::Array{ComplexF64,2}, C::Array{ComplexF64,2})
    A = overlap_gamma(Ctilde, Ctilde)
    B = overlap_gamma(C, Ctilde)
    X = 0.5*(I - A)
    Xnew = similar(X)
    for iter in 1:100
        Xnew[:,:] = 0.5*( I - A + X*(I - B) + (I-B)*X - X*X )
        #println("real X - Xnew")
        #display(real(X - Xnew)); println()
        #println("imag X-Xnew")
        #display(imag(X - Xnew)); println()
        ΔX = norm(X-Xnew)
        println("norm X-Xnew = ", ΔX)
        if ΔX < 1e-10
            println("Find X: converged in iter: ", iter)
            return Xnew
        end
        @views X[:,:] = Xnew[:,:]
    end
    println("ERROR: X is not converged")
    exit()
    return X
end


function main( init_func; fnametrj="TRAJ.xyz", fnameetot="ETOT.dat" )

    dt_fs = 0.125
    # Time step, in Ha atomic unit
    dt = dt_fs*10e-16/AU_SEC
    println("dt (au) = ", dt)
    println("dt (ps) = ", dt*AU_PS)

    Ham = init_func()

    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates
    atoms = Ham.atoms
    Natoms = atoms.Natoms
    m = atoms.masses
    atm2species = atoms.atm2species
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    println(Ham.atoms.masses)

    # Ionic variables
    p = zeros(Float64,3,Natoms)
    dr = zeros(Float64,3,Natoms)
    v = zeros(Float64,3,Natoms)
    vtilde = zeros(Float64,3,Natoms)

    # Electronic variables (other than psis)
    μ = 400.0 # fictitious mass of electron

    filetraj = open(fnametrj, "w")
    fileetot = open(fnameetot, "w")

    C = randn_BlochWavefuncGamma(Ham)
    energies, forces = minimize_electrons!(Ham, C, etot_conv_thr=1e-8)

    Etot = sum(energies)
    Ekin_elec = 0.0
    Ekin_ions = 0.0  # assume initial velocities are zeroes
    Etot_conserved = Etot + Ekin_ions + Ekin_elec

    # 0-th iteration, initial condition
    iter = 0
    write_traj_files(
        filetraj, fileetot,
        atoms, forces,
        Etot_conserved, Etot, Ekin_ions, Ekin_elec,
        dt, iter
    )


    @assert Ham.electrons.Nspin == 1
    F_elec = zeros(ComplexF64,size(C.data[1]))
    Ctilde = zeros(ComplexF64,size(C.data[1]))
    dCtilde = zeros(ComplexF64,size(C.data[1]))
    dC = zeros(ComplexF64,size(C.data[1]))
    dCpr = zeros(ComplexF64,size(C.data[1]))

    Rhoe = similar(Ham.rhoe)

    calc_grad!(Ham, C.data[1], F_elec)
    F_elec = -F_elec

    NiterMax = 5000
    for iter = 1:NiterMax

        # Update atomic positions
        for ia in 1:Natoms
            isp = atm2species[ia]
            for i in 1:3
                vtilde[i,ia] = v[i,ia] + 0.5*dt*forces[i,ia]/m[isp]
                dr[i,ia] = dt*vtilde[i,ia]
            end
        end
        update_positions!( Ham, Ham.atoms.positions + dr )

        #
        dCtilde[:] = dC + 0.5*dt*F_elec/μ
        Ctilde[:] = C.data[1] + dt*dCtilde
        X = calc_X_matrix( Ctilde, C.data[1] )
        C.data[1][:,:] = Ctilde + C.data[1]*X
        
        println("Test ortho: ")
        println(dot(C,C))
        println( dot_gamma(C.data[1][:,1], C.data[1][:,1]) )
        println( dot_gamma(C.data[1][:,1], C.data[1][:,2]) )
        println( dot_gamma(C.data[1][:,1], C.data[1][:,3]) )
        #exit()

        Rhoe[:,:] = calc_rhoe(Ham, C)
        update!(Ham, Rhoe)
        println("integ Rhoe = ", sum(Rhoe)*dVol)

        energies = calc_energies(Ham, C)
        Etot = sum(energies)

        forces[:,:] = calc_forces(Ham, C)
        calc_grad!(Ham, C.data[1], F_elec)
        F_elec = -F_elec

        # Update ions' velocities
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

        # Update electrons' velocities
        dCpr[:] = Ctilde + 0.5*dt*F_elec/μ
        Q = overlap_gamma(C.data[1], dCpr)
        Y = -0.5*(Q + Q')
        dC[:] = dCpr + C.data[1]*Y

        # Orthogonality condition on orbital velocities
        println("test1 = ", dot_gamma(dC, C.data[1]))
        println("test2 = ", dot_gamma(C.data[1], dC))

        Ekin_elec = 0.5*μ*real(dot_gamma(dC,dC))
        println("Ekin_elec = ", Ekin_elec)

        Etot_conserved = Etot + Ekin_ions + Ekin_elec

        write_traj_files(
            filetraj, fileetot,
            atoms, forces,
            Etot_conserved, Etot, Ekin_ions, Ekin_elec,
            dt, iter
        )


        if Ekin_elec > 1e-3
            println("Quenching ...")
            energies, forces[:] = minimize_electrons!(Ham, C)
            dC[:] .= 0.0
            calc_grad!(Ham, C.data[1], F_elec)
            F_elec[:] = -F_elec
        end


        #@printf("\nMD Iter = %3d, Etot = %18.10f\n", iter, sum(energies))
        #println("Forces = ")
        #for ia in 1:Natoms
        #    isp = atm2species[ia]
        #    atsymb = Ham.atoms.SpeciesSymbols[isp]
        #    @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
        #        forces[1,ia], forces[2,ia], forces[3,ia])
        #end

    end

    close(filetraj)
    close(fileetot)

end

#main(init_Ham_H2O, fnametrj="TRAJ_H2O_v4.xyz", fnameetot="ETOT_H2O_v4.dat")

main(init_Ham_CO2,
    fnametrj="TRAJ_CO2_cpmd.xyz",
    fnameetot="ETOT_CO2_cpmd.dat"
)

#main(init_Ham_Si8,
#    fnametrj="TRAJ_Si8_cpmd.xyz",
#    fnameetot="ETOT_Si8_cpmd.dat"
#)

#main(init_Ham_CO2,
#    fnametrj="TRAJ_CO2_step10_extrap2nd.xyz",
#    fnameetot="ETOT_CO2_step10_extrap2nd.dat"
#)