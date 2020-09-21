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

function initial_hessian( Natoms )
    #h = 70/(2*Ry2eV) * (ANG2BOHR^2)
    # Convert to Ha/(bohr^2)
    h = 70*eV2Ha/ANG2BOHR^2
    H = diagm( 0 => h*ones(3*Natoms) )
    return H
end
    
function update_hessian( H_old, r, f, r0, f0 )
    Natoms = size(r,2)
    dr = r - r0
    # FIXME: need this?
    if maximum( abs.(dr) ) < 1e-7
        return diagm( 0 => 70*ones(3*Natoms) )
    end
    df = f - f0
    a = dot(dr, df)
    dg = H_old*dr
    b = dot(dr, dg)
    return H_old - (df*df')/a - (dg*dg')/b
end

function run_pwdft_jl!( Ham, psis )
    KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=200 )
    forces = calc_forces( Ham, psis )
    return sum(Ham.energies), forces
end

function do_bfgs( molname::String; ecutwfc=15.0 )
    
    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", molname*".xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc)

    Natoms = Ham.atoms.Natoms

    psis = randn_BlochWavefuncGamma(Ham)
    energies, forces = run_pwdft_jl!(Ham, psis)

    filetraj = open("TRAJ_"*molname*".xyz", "w")

    @printf(filetraj, "%d\n\n", Natoms)
    atsymbs = atoms.atsymbs
    τ = atoms.positions*BOHR2ANG # convert to Angstrom
    FF = forces*(1/BOHR2ANG) # convert to Ha/Angstrom
    for ia in 1:Natoms
        @printf(filetraj, "%s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
            atsymbs[ia], τ[1,ia], τ[2,ia], τ[3,ia],
            forces[1,ia], forces[2,ia], forces[3,ia])
    end
    
    MAXSTEP = 0.04*ANG2BOHR

    f = vec(copy(forces))
    r = vec(copy(Ham.atoms.positions))
    r0 = zeros(size(r))
    f0 = zeros(size(f))

    H = initial_hessian(Natoms)

    NiterMax = 15
    for iter = 1:NiterMax

        println("Hessian = ")
        display(H); println()

        omega, V = eigen( Symmetric(H) )
        dr = V * (V'*f ./ abs.(omega))
        steplengths = sqrt.(sum( dr.^2, dims=1 ))
        maxsteplength = maximum(steplengths)
        if maxsteplength >= MAXSTEP
            println("Scaling dr")
            dr = dr * MAXSTEP / maxsteplength
        end

        r0 = copy(r)
        f0 = copy(f)
        H_old = copy(H)

        r[:] = r[:] + dr[:]
        #Ham.atoms.positions = copy(reshape(r,(3,Natoms)))
        update_positions!( Ham, reshape(r,3,Natoms) )

        energies_old = energies
        energies, forces = run_pwdft_jl!(Ham, psis)

        @printf(filetraj, "%d\n\n", Natoms)
        atsymbs = atoms.atsymbs
        τ = atoms.positions*BOHR2ANG # convert to Angstrom
        FF = forces*(1/BOHR2ANG) # convert to Ha/Angstrom
        for ia in 1:Natoms
            @printf(filetraj, "%s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                atsymbs[ia], τ[1,ia], τ[2,ia], τ[3,ia],
                forces[1,ia], forces[2,ia], forces[3,ia])
        end

        fmax = norm(forces)

        @printf("\nIter = %3d, Etot = %18.10f, fmax = %10.5e\n", iter, sum(energies), fmax)
        println("Forces = ")
        display(forces'); println()
        println("dr = ")
        display(reshape(dr,(3,Natoms))'); println()
        println("r  =")
        display(Ham.atoms.positions'); println()
        println("r (in Angstrom) =")
        display(Ham.atoms.positions' / ANG2BOHR); println()

        if fmax < 0.01*eV2Ha/ANG2BOHR
            println("BFGS is converged")
            break
        end

        f = vec(copy(forces))
        H = update_hessian( H_old, r, f, r0, f0 )

    end

    close(filetraj)

end

function main()
    Nargs = length(ARGS)
    if Nargs >= 1
        molname = ARGS[1]
    else
        molname = "H2O"
    end
    do_bfgs(molname)
end

main()
