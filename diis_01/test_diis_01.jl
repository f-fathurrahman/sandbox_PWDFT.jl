using Printf
using LinearAlgebra
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("create_Ham.jl")

# Similar to update!, but without copying Rhoe
function update_pots_from_rhoe!( Ham::Hamiltonian )

    @assert Ham.rhoe_core == nothing
    @assert Ham.electrons.Nspin == 1

    # Ham.rhoe is assumed to be updated elsewhere
    Rhoe = Ham.rhoe

    pw = Ham.pw
    xc_calc = Ham.xc_calc
    V_Hartree = Ham.potentials.Hartree
    V_XC = Ham.potentials.XC
    V_Total = Ham.potentials.Total
    V_Ps_loc = Ham.potentials.Ps_loc

    Rhoe_tot = dropdims(Rhoe, dims=2) # should be sum if Nspin != 1
    Poisson_solve!(pw, Rhoe_tot, V_Hartree)

    if Ham.xcfunc == "PBE"
        @views V_XC[:,1] .= calc_Vxc_PBE( xc_calc, pw, Rhoe )
    else
        # use VWN
        @views V_XC[:,1] .= calc_Vxc_VWN( xc_calc, Rhoe )
    end
    
    # Update total potentials
    Npoints = size(Rhoe,1)
    for ip in 1:Npoints
        V_Total[ip,1] = V_Ps_loc[ip] + V_Hartree[ip] + V_XC[ip,1]  
    end

    return
end

function main()

    #Ham = create_Ham_CO()
    #Ham = create_Ham_NH3()
    #Ham = create_Ham_N2()
    Ham = create_Ham_G2_mols(molname="H2O")

    #KS_solve_Emin_PCG!( Ham )
    #return

    @assert Ham.pw.gvecw.kpoints.Nkpt == 1

    Npoints = prod(Ham.pw.Ns)
    CellVolume = Ham.pw.CellVolume
    dVol = CellVolume/Npoints
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin

    Rhoe = Ham.rhoe

    Rhoe[:,1] = guess_rhoe( Ham )

    println("sum Rhoe = ", sum(Ham.rhoe)*dVol)

    update_pots_from_rhoe!(Ham)

    evals = Ham.electrons.ebands
    psiks = rand_BlochWavefunc(Ham)

    # Initial diagonalization to get psiks
    evals[:,:] = diag_LOBPCG!( Ham, psiks, NiterMax=10 )


    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # Calculate energy at this psi
    energies = calc_energies(Ham, psiks)
    Ham.energies = energies
    Etot = sum(energies)

    Etot_old = 0.0  # declare the variable

    Ndiis = 5
    psiks_diis = Vector{BlochWavefunc}(undef, Ndiis)
    for i in 1:Ndiis
        psiks_diis[i] = zeros_BlochWavefunc(Ham)
    end
    #
    for ikspin in 1:length(psiks)
        psiks_diis[1][ikspin][:,:] .= psiks[ikspin][:,:]
    end

    g_diis = Vector{BlochWavefunc}(undef, Ndiis)
    Kg_diis = Vector{BlochWavefunc}(undef, Ndiis)
    for i in 1:Ndiis
        g_diis[i] = zeros_BlochWavefunc(Ham)
        Kg_diis[i] = zeros_BlochWavefunc(Ham)
    end

    _calc_all_grads!(Ham, psiks_diis[1], g_diis[1], Kg_diis[1])
    for i in 2:Ndiis
        Etot_old = Etot
        println("\nSteepest descent step = ", i)
        Etot = _steepest_descent_step!(
            Ham, psiks_diis[i-1],
            g_diis[i-1], Kg_diis[i-1],
            psiks_diis[i]
        )
        @printf("iterSD %5d %18.10f %18.10e\n", i, Etot, abs(Etot - Etot_old))
        # Calculate new gradients
        _calc_all_grads!(Ham, psiks_diis[i], g_diis[i], Kg_diis[i])
        #
    end

    d_diis = _calc_d_diis(Kg_diis)
    println("real d_diis = ", real(d_diis))
    println("imag d_diis = ", imag(d_diis))
    println("sum d_diis = ", sum(d_diis))


    # Run more SD steps
    NMoreSDSteps = 0
    for iterSD in Ndiis+1:NMoreSDSteps

        println("\nBegin iterSD = ", iterSD)

        # FIXME: Shift, until some convergence is reached then start using DIIS
        for i in 2:Ndiis
            psiks_diis[i-1][1][:,:] = psiks_diis[i][1][:,:]
            g_diis[i-1][1][:,:] = g_diis[i][1][:,:]
            Kg_diis[i-1][1][:,:] = Kg_diis[i][1][:,:]
        end
        # *_diis[1][1][:,:] are discarded
        #
        #
        Etot_old = Etot
        Etot = _steepest_descent_step!(
            Ham, psiks_diis[Ndiis-1],
            g_diis[Ndiis-1], Kg_diis[Ndiis-1],
            psiks_diis[Ndiis]
        )
        # Calculate Kg_diis for idiis = Ndiis
        _calc_all_grads!(Ham, psiks_diis[Ndiis], g_diis[Ndiis], Kg_diis[Ndiis])

        # Check DIIS
        d_diis = _calc_d_diis(Kg_diis)
        println("real d_diis = ", real(d_diis))
        println("imag d_diis = ", imag(d_diis))
        println("sum d_diis = ", sum(d_diis))

        @printf("iterSD %5d %18.10f %18.10e\n", iterSD, Etot, abs(Etot - Etot_old))

    end


    MAX_ITER_DIIS = 100 # set to zero to disable DIIS

    # New psiks and g, and Kg
    ik = 1 # limit to Nkspin = 1
    ispin = 1
    ikspin = 1
    #
    g_new = zeros_BlochWavefunc(Ham)
    Kg_new = zeros_BlochWavefunc(Ham)
    psiks_new = zeros_BlochWavefunc(Ham)

    for iterDIIS in 1:MAX_ITER_DIIS

        # DIIS procedure starts here
        d_diis = _calc_d_diis(Kg_diis)

        println("real d_diis = ", real(d_diis))
        println("imag d_diis = ", imag(d_diis))
        println("sum d_diis = ", sum(d_diis))

        # Approximate g
        fill!(g_new[ikspin], 0.0)
        for i in 1:Ndiis
            g_new[ikspin][:,:] += d_diis[i]*g_diis[i][ikspin][:,:]
        end
        Kg_new[ikspin][:,:] .= Kprec(ik, Ham.pw, g_new[ikspin])

        # New psiks
        fill!(psiks_new[ikspin], 0.0)
        for i in 1:Ndiis
            psiks_new[ikspin][:,:] += d_diis[i]*psiks_diis[i][ikspin][:,:]
        end
        psiks_new[ikspin][:,:] -= Kg_new[ikspin][:,:] # Note that the sign is different

        # Don't forget to orthogonalize
        ortho_sqrt!( psiks_new[ikspin] )

        # New Hamiltonian, calculate energies
        calc_rhoe!(Ham, psiks_new, Rhoe)
        println("DIIS: Integ Rhoe = ", sum(Rhoe)*dVol)
        update_pots_from_rhoe!(Ham)

        Ham.energies = calc_energies( Ham, psiks_new )
        Etot_old = Etot
        Etot = sum(Ham.energies)
        diffE = abs(Etot - Etot_old)
    
        @printf("iterDIIS %5d %18.10f %18.10e\n", iterDIIS, Etot, diffE)

        if diffE <= 1e-6
            println("DIIS converged")
            break
        end

        # Calculate actual g using psiks_new
        for i in 1:Ndiis
            g_new[ikspin][:,:] .= calc_grad(Ham, psiks_new[ikspin])
        end
        Kg_new[ikspin][:,:] .= Kprec(ik, Ham.pw, g_new[ikspin])
        # FIXME: use calc_all_grads

        for i in 2:Ndiis
            psiks_diis[i-1][ikspin][:,:] = psiks_diis[i][ikspin][:,:]
            g_diis[i-1][ikspin][:,:] = g_diis[i][ikspin][:,:]
            Kg_diis[i-1][ikspin][:,:] = Kg_diis[i][ikspin][:,:]
        end
        psiks_diis[Ndiis][ikspin][:,:] = psiks_new[ikspin][:,:]
        g_diis[Ndiis][ikspin][:,:] = g_new[ikspin][:,:] # use the actual g_new?
        Kg_diis[Ndiis][ikspin][:,:] = Kg_new[ikspin][:,:] # gnew ?

    end

    return

end


function _calc_d_diis(Kg_diis)
    
    Ndiis = length(Kg_diis)
    Nkspin = length(Kg_diis[1])

    @assert Nkspin == 1

    # FIXME: B matrix and d-vector have Nkspin dimension
    B = ones(ComplexF64, Ndiis+1, Ndiis+1)
    B[Ndiis+1,Ndiis+1] = 0

    ikspin = 1
    for k in 1:Ndiis, l in 1:Ndiis
        B[k,l] = dot(Kg_diis[k][ikspin], Kg_diis[l][ikspin])
    end

    println("\nDIIS matrix B:\n")
    display(real(B)); println() # only show real parts for simplicity
    #display(imag(B)); println()

    rhs_vec = zeros(ComplexF64, Ndiis+1)
    rhs_vec[Ndiis+1] = 1

    d_full = B\rhs_vec

    return d_full[1:Ndiis]
end

function _calc_all_grads!(Ham, psiks, g, Kg)
    #
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    #
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        g[i][:,:] .= calc_grad( Ham, psiks[i] )
        Kg[i][:,:] .= Kprec( Ham.ik, Ham.pw, g[i] )
    end
    return
end

# g and Kg is calculated before calling this function
function _steepest_descent_step!(Ham, psiks, g, Kg, psiks_new)
    
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin

    Rhoe = Ham.rhoe # alias

    d    = zeros_BlochWavefunc(Ham)
    psic = zeros_BlochWavefunc(Ham)
    gt   = zeros_BlochWavefunc(Ham)
    
    α_t = 3e-5

    β = zeros(Nkspin)
    α = zeros(Nkspin)

    for ispin in 1:Nspin, ik in 1:Nkpt
        i = ik + (ispin - 1)*Nkpt
        d[i][:,:] = -Kg[i][:,:]
        psic[i] = ortho_sqrt(psiks[i] + α_t*d[i])
    end
        
    calc_rhoe!( Ham, psic, Rhoe )
    update_pots_from_rhoe!(Ham)

    # Update using line minimization
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        #
        gt[i] = calc_grad(Ham, psic[i])
        #
        denum = real(sum(conj(g[i]-gt[i]).*d[i]))
        if denum != 0.0
            α[i] = abs( α_t*real(sum(conj(g[i]).*d[i]))/denum )
        else
            α[i] = 0.0
        end
        # Update wavefunction
        psiks_new[i][:,:] = psiks[i] + α[i]*d[i]
        ortho_sqrt!(psiks_new[i])
    end

    # Update Rhoe and potentials
    calc_rhoe!(Ham, psiks_new, Rhoe)

    update_pots_from_rhoe!(Ham)

    Ham.energies = calc_energies( Ham, psiks_new )
    Etot = sum(Ham.energies)

    return Etot
end

main()