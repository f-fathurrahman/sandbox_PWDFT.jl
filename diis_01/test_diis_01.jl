using Printf
using LinearAlgebra
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

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

    Ham = create_Ham_CO()

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

    for i in 2:Ndiis
        println("\nSteepest descent step = ", i)
        _steepest_descent_step!(
            Ham, psiks_diis[i-1],
            g_diis[i-1], Kg_diis[i-1],
            psiks_diis[i]
        )
    end

    # Calculate Kg_diis for idiis = Ndiis
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        g_diis[Ndiis][i][:,:] .= calc_grad( Ham, psiks[i] )
        Kg_diis[Ndiis][i][:,:] .= Kprec( Ham.ik, Ham.pw, g_diis[Ndiis][i] )
    end

    for NNNN in 1:2
        # FIXME: Shift, until some convergence is reached then start using DIIS
        for i in 2:Ndiis
            psiks_diis[i-1][1][:,:] = psiks_diis[i][1][:,:]
            g_diis[i-1][1][:,:] = g_diis[i][1][:,:]
            Kg_diis[i-1][1][:,:] = Kg_diis[i][1][:,:]
        end
        #
        _steepest_descent_step!(
            Ham, psiks_diis[Ndiis-1],
            g_diis[Ndiis-1], Kg_diis[Ndiis-1],
            psiks_diis[Ndiis]
        )
        # Calculate Kg_diis for idiis = Ndiis
        for ispin in 1:Nspin, ik in 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            i = ik + (ispin - 1)*Nkpt
            g_diis[Ndiis][i][:,:] .= calc_grad( Ham, psiks[i] )
            Kg_diis[Ndiis][i][:,:] .= Kprec( Ham.ik, Ham.pw, g_diis[Ndiis][i] )
        end
    end


    d_diis = _calc_d_diis(Kg_diis)

    println("sum d_diis = ", sum(d_diis))

    # New psiks and g, and Kg
    ik = 1 # limit to Nkspin = 1
    ispin = 1
    ikspin = 1
    #
    g_new = zeros_BlochWavefunc(Ham)
    Kg_new = zeros_BlochWavefunc(Ham)
    # Approximate g
    for i in 1:Ndiis
        g_new[ikspin][:,:] += d_diis[i]*g_diis[i][ikspin][:,:]
    end
    Kg_new[ikspin][:,:] .= Kprec(ik, Ham.pw, g_new[ikspin])

    psiks_new = zeros_BlochWavefunc(Ham)
    for i in 1:Ndiis
        psiks_new[ikspin][:,:] += d_diis[i]*psiks_diis[i][ikspin][:,:]
    end
    psiks_new[ikspin][:,:] += Kg_new[ikspin][:,:]

    # Don't forget to orthogonalize
    ortho_sqrt!( psiks_new[ikspin] )

    # New Hamiltonian, calculate energies
    calc_rhoe!(Ham, psiks_new, Rhoe)
    update_pots_from_rhoe!(Ham)

    Ham.energies = calc_energies( Ham, psiks_new )
    Etot = sum(Ham.energies)
    println("Etot = ", Etot)


    println("Pass here")

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
    display(real(B)); println()
    display(imag(B)); println()

    rhs_vec = zeros(ComplexF64, Ndiis+1)
    rhs_vec[Ndiis+1] = 1

    d_full = B\rhs_vec
    println(d_full[1:Ndiis])

    return d_full[1:Ndiis]
end


function _steepest_descent_step!(Ham, psiks, g, Kg, psiks_new)
    
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Rhoe = Ham.rhoe

    d    = zeros_BlochWavefunc(Ham)
    psic = deepcopy(g)
    gt   = deepcopy(g)
    
    α_t = 3e-5

    β = zeros(Nkspin)
    α = zeros(Nkspin)

    for ispin in 1:Nspin, ik in 1:Nkpt

        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt

        g[i][:,:] .= calc_grad( Ham, psiks[i] )
        Kg[i][:,:] .= Kprec( Ham.ik, Ham.pw, g[i] )

        d[i] = -Kg[i]

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
    println("Etot = ", Etot)

    return
end

main()