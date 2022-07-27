using Printf
using Random
using LinearAlgebra
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("../get_default_psp.jl")
include("update_positions.jl")

Random.seed!(1234)

include("additional_constants.jl")
include("init_Ham.jl")
include("utils_md.jl")

# Initial minimization
function minimize_electrons!( Ham, psiks; NiterMax=200, etot_conv_thr=1e-8 )
    KS_solve_Emin_PCG!( Ham, psiks,
        skip_initial_diag=true, etot_conv_thr=etot_conv_thr,
        NiterMax=NiterMax
    )
    forces = calc_forces( Ham, psiks )
    return sum(Ham.energies), forces
end

# μ is fictitious electron mass
function calc_elec_energies!(Ham, psiks, μ, dpsiks)
    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, psiks)
    energies = calc_energies(Ham, psiks)
    Etot = sum(energies)
    Ekin_elec = 0.5*μ*dot(dpsiks,dpsiks)
    return Etot, Ekin_elec
end 

# RATTLE algorithm
function calc_X_matrix(Ctilde::Array{ComplexF64,2}, C::Array{ComplexF64,2})
    A = Ctilde' * Ctilde
    B = C' * Ctilde
    X = 0.5*(I - A)
    Xnew = similar(X)
    for iter in 1:100
        Xnew[:,:] = 0.5*( I - A + X*(I - B) + (I - B)*X - X*X )
        ΔX = norm(X-Xnew)
        println("norm X-Xnew = ", ΔX)
        if ΔX < 1e-10
            println("Find X: converged in iter: ", iter)
            return Xnew
        end
        @views X[:,:] = Xnew[:,:]
    end
    error("ERROR: X is not converged")
    return X
end

# Need this?
struct CarParrinelloMDParameters
    dt_fs::Float64
    dt::Float64
    μ::Float64
end

function main( init_func; fnametrj="TRAJ.xyz", fnameetot="ETOT.dat" )

    dt_fs = 0.125
    #dt_fs = 0.3
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

    # Electronic variables (other than psiks)
    #μ = 400.0 # fictitious mass of electron
    μ = 800.0 # fictitious mass of electron

    filetraj = open(fnametrj, "w")
    fileetot = open(fnameetot, "w")

    psiks = rand_BlochWavefunc(Ham)
    energies, forces = minimize_electrons!(Ham, psiks, etot_conv_thr=1e-8)

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
    @assert Ham.pw.gvecw.kpoints.Nkpt == 1

    C = psiks[1]

    F_elec = zeros(ComplexF64,size(C))
    Ctilde = zeros(ComplexF64,size(C))
    dCtilde = zeros(ComplexF64,size(C))
    dC = zeros(ComplexF64,size(C))
    dCpr = zeros(ComplexF64,size(C))

    Rhoe = similar(Ham.rhoe)

    #calc_grad!(Ham, C, F_elec)
    F_elec = -calc_grad(Ham, C)
    #F_elec = -F_elec

    NiterMax = 1

    #
    for iter in 1:NiterMax

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
        Ctilde[:] = C + dt*dCtilde
        X = calc_X_matrix( Ctilde, C )
        C[:,:] = Ctilde + C*X

        println("Test ortho: ")
        println(dot(C,C))
        println( dot(C[:,1], C[:,1]) )
        println( dot(C[:,1], C[:,2]) )
        println( dot(C[:,1], C[:,3]) )
        #exit()

        Rhoe[:,:] = calc_rhoe(Ham, psiks)
        update!(Ham, Rhoe)
        println("integ Rhoe = ", sum(Rhoe)*dVol)

        energies = calc_energies(Ham, psiks)
        Etot = sum(energies)

        forces[:,:] = calc_forces(Ham, psiks)
        F_elec = -calc_grad(Ham, C)
        #calc_grad!(Ham, C, F_elec)
        #F_elec = -F_elec

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
        Q = C' * dCpr
        Y = -0.5*(Q + Q') # Marx Eq. 3.135
        dC[:] = dCpr + C*Y
        # Y can be ontained without iteration
        # The velocity constraint condition is satisfied exacatly,
        # up to machine precision, at teach time step
        # Marx Sec. 3.7.1

        # Orthogonality condition on orbital velocities
        println("test1 dot dC,C = ", dot(dC, C))
        println("test2 dot C,dC = ", dot(C, dC))

        Ekin_elec = 0.5*μ*real(dot(dC,dC))
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
            energies, forces[:] = minimize_electrons!(Ham, psiks)
            dC[:] .= 0.0
            #calc_grad!(Ham, C, F_elec)
            F_elec[:,:] = calc_grad(Ham, C)
            #F_elec[:] = -F_elec
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

#main(init_Ham_CO2,
#    fnametrj="TRAJ_CO2_cpmd.xyz",
#    fnameetot="ETOT_CO2_cpmd.dat"
#)

main(init_Ham_Si8,
    fnametrj="TRAJ_Si8_cpmd.xyz",
    fnameetot="ETOT_Si8_cpmd.dat"
)

#main(init_Ham_CO2,
#    fnametrj="TRAJ_CO2_step10_extrap2nd.xyz",
#    fnameetot="ETOT_CO2_step10_extrap2nd.dat"
#)