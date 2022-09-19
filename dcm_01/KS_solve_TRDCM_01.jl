struct DCMBlockSets
    set1::UnitRange{Int64}
    set2::UnitRange{Int64}
    set3::UnitRange{Int64}
    set4::UnitRange{Int64}
    set5::UnitRange{Int64}
end

function DCMBlockSets(Nstates)
    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates
    set3 = 2*Nstates+1:3*Nstates
    set4 = Nstates+1:3*Nstates
    set5 = 1:2*Nstates
    return DCMBlockSets(set1, set2, set3, set4, set5)
end



"""
Solves Kohn-Sham problem using trust-region direct constrained minimization
(DCM) as described by Prof. Chao Yang.
"""
function KS_solve_TRDCM_01!(
    Ham::Hamiltonian;
    NiterMax=100,
    startingwfc=:random,
    verbose=true,
    print_final_ebands=false,
    print_final_energies=true,
    etot_conv_thr=1e-6
)


	pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    ΔV = CellVolume/Npoints
    electrons = Ham.electrons
    Focc = electrons.Focc
    Nstates = electrons.Nstates
    Nocc = electrons.Nstates_occ
    Nelectrons = electrons.Nelectrons

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = electrons.Nspin

    Nkspin = Nkpt*Nspin

    psiks = BlochWavefunc(undef,Nkspin)

    #
    # Initial wave function
    #
    if startingwfc == :random
        # generate random BlochWavefunc
        psiks = rand_BlochWavefunc( Ham )
    end

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)
    #if startingrhoe == :gaussian
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
    #end
    #calc_rhoe!( Ham, psiks, Rhoe )
    update!(Ham, Rhoe)

    evals = zeros(Float64,Nstates,Nkspin)
    
    # Starting eigenvalues and psi
    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        evals[:,i] = diag_LOBPCG!( Ham, psiks[i], verbose_last=false, NiterMax=10 )
    end

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # XXX: Update Rhoe?
    #calc_rhoe!( Ham, psiks, Rhoe )
    #update!(Ham, Rhoe)

    #
    Ham.energies = calc_energies( Ham, psiks )
    Etot = sum(Ham.energies)
    Etot_old = Etot

    # subspace
    Y = BlochWavefunc(undef,Nkspin)
    R = BlochWavefunc(undef,Nkspin)
    P = BlochWavefunc(undef,Nkspin)
    #
    G = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    T = Array{Array{ComplexF64,2},1}(undef,Nkspin) # XXX: can be made real?
    B = Array{Array{ComplexF64,2},1}(undef,Nkspin) # XXX: can be made real?
    A = Array{Array{ComplexF64,2},1}(undef,Nkspin) # XXX: can be made real?
    C = Array{Array{ComplexF64,2},1}(undef,Nkspin) # XXX: should be complex
    for ispin in 1:Nspin, ik in 1:Nkpt
        i = ik + (ispin - 1)*Nkpt
        Y[i] = zeros( ComplexF64, Ngw[ik], 3*Nstates )
        R[i] = zeros( ComplexF64, Ngw[ik], Nstates )
        P[i] = zeros( ComplexF64, Ngw[ik], Nstates )
        G[i] = zeros( ComplexF64, 3*Nstates, 3*Nstates ) # also used to hold eigenvector of reduced problem
        T[i] = zeros( ComplexF64, 3*Nstates, 3*Nstates )
        B[i] = zeros( ComplexF64, 3*Nstates, 3*Nstates )
        A[i] = zeros( ComplexF64, 3*Nstates, 3*Nstates )
        C[i] = zeros( ComplexF64, 3*Nstates, 3*Nstates )
    end

    # array for saving eigenvalues of subspace problem
    D = zeros(Float64,3*Nstates,Nkspin)

    #XXX use plain 3d-array for G, T, and B ?

    dcm_block_sets = DCMBlockSets(Nstates)
    set1 = dcm_block_sets.set1
    set2 = dcm_block_sets.set2
    set3 = dcm_block_sets.set3
    set4 = dcm_block_sets.set4
    set5 = dcm_block_sets.set5

    MaxInnerSCF = 3
    MAXTRY = 10
    FUDGE = 1e-12
    SMALL = 1e-12

    σ = zeros(Float64,Nkspin)
    gapmax = zeros(Float64,Nkspin)

    diffE = 1.0 # for evaluating convergence

    Nconverges = 0

    for iter in 1:NiterMax
        
        _dcm_construct_Y_and_G!(Ham, iter, dcm_block_sets, psiks, R, T, B, P, Y, G)
        
        @printf("TRDCM iter: %3d\n", iter)

        σ[:] .= 0.0  # reset σ to zero at the beginning of inner SCF iteration
        numtry = 0
        
        Etot_innerscf = sum(Ham.energies)
        Etot_innerscf_old = Etot_innerscf

        println("Etot_innerscf = ", Etot_innerscf)

        for iterscf in 1:MaxInnerSCF
            
            # psiks is also updated            
            _dcm_project_nonlinear_pot!(Ham, iter, dcm_block_sets, psiks, σ, evals, A, B, C, D, T, Y, G)

            calc_rhoe!( Ham, psiks, Rhoe )
            #println("integ Rhoe = ", sum(Rhoe)*ΔV)
            update!( Ham, Rhoe )

            # Calculate energies once again
            Ham.energies = calc_energies( Ham, psiks )
            Etot_innerscf = sum(Ham.energies)

            println("Etot_innerscf = ", Etot_innerscf)

            if Etot_innerscf > Etot_innerscf_old

                @printf("TRDCM: %f > %f: Trust region will be imposed\n", Etot_innerscf, Etot_innerscf_old)

                # Total energy is increased, impose trust region
                # Do this for all kspin

                for i in 1:Nkspin

                    if iter == 1
                        gaps = D[2:2*Nstates,i] - D[1:2*Nstates-1,i]
                        gapmax[i] = maximum(gaps)
                    else
                        gaps = D[2:3*Nstates] - D[1:3*Nstates-1]
                        gapmax[i] = maximum(gaps)
                    end
                    gap0 = D[Nocc+1,i] - D[Nocc,i]

                    while (gap0 < 0.9*gapmax[i]) && (numtry < MAXTRY)
                        println("Increase σ to fix gap0: numtry = ", numtry)
                        @printf("gap0 : %10.5f < %10.5f\n", gap0, 0.9*gapmax[i])
                        if abs(σ[i]) < SMALL # approx for σ == 0.0
                            # initial value for σ
                            σ[i] = 2*gapmax[i]
                        else
                            σ[i] = 2*σ[i]
                        end
                        @printf("fix gap0: i = %3d, σ = %18.10f\n", i, σ[i])
                        #
                        if iter > 1
                            D[:,i], G[i] = eigen( A[i] - σ[i]*C[i], B[i] )
                            gaps = D[2:2*Nstates,i] - D[1:2*Nstates-1,i]
                        else
                            D[set5,i], G[i][set5,set5] =
                            eigen(A[i][set5,set5] - σ[i]*C[i][set5,set5], B[i][set5,set5] )
                            gaps = D[2:3*Nstates,i] - D[1:3*Nstates-1,i]
                        end
                        gapmax[i] = maximum(gaps)
                        gap0 = D[Nocc+1,i] - D[Nocc,i]
                        numtry = numtry + 1
                    end
                end # Nkspin

            end # if Etot > Etot0

            #println("σ = ", σ)
            numtry = 0  # reset numtry for this while loop

            while (Etot_innerscf > Etot_innerscf_old) &
                  #(abs(Etot-Etot0) > FUDGE*abs(Etot0)) &
                  (abs(Etot_innerscf-Etot_innerscf_old) > etot_conv_thr*0.01) &
                  (numtry < MAXTRY)

                @printf("Increase σ part 2: %18.10f > %18.10f\n", Etot_innerscf, Etot_innerscf_old)
                @printf("Increase σ part 2: diff = %18.10e\n", Etot_innerscf - Etot_innerscf_old)
                #
                # update wavefunction
                #
                for i = 1:Nkspin
                    if iter > 1
                        psiks[i] = Y[i]*G[i][:,set1]
                        #ortho_sqrt!(psiks[i])
                    else
                        psiks[i] = Y[i][:,set5]*G[i][set5,set1]
                        #ortho_sqrt!(psiks[i])
                    end
                end

                calc_rhoe!( Ham, psiks, Rhoe )
                update!( Ham, Rhoe )
            
                # Calculate energies once again
                Ham.energies = calc_energies( Ham, psiks )
                Etot_innerscf = sum(Ham.energies)
                
                if Etot_innerscf > Etot_innerscf_old
                    
                    println("Increase σ part 2")
                    
                    for i in 1:Nkspin
                        if abs(σ[i]) > SMALL # σ is not 0
                            σ[i] = 2*σ[i]
                        else
                            σ[i] = 1.2*gapmax[i]
                        end
                        @printf("ikspin = %d σ = %f\n", i, σ[i])
                        if iter > 1
                            D[:,i], G[i] = eigen(A[i] - σ[i]*C[i], B[i])
                        else
                            D[set5,i], G[i][set5,set5] =
                            eigen(A[i][set5,set5] - σ[i]*C[i][set5,set5], B[i][set5,set5])
                        end
                    end

                end
                numtry = numtry + 1  # outside i loop
            end # while

            Etot_innerscf_old = Etot_innerscf
            
        end # end of inner SCF iteration


        for i in 1:Nkspin
            if iter > 1
                psiks[i] = Y[i]*G[i][:,set1]
                #ortho_sqrt!(psiks[i])
            else
                psiks[i] = Y[i][:,set5]*G[i][set5,set1]
                #ortho_sqrt!(psiks[i])
            end
        end
        calc_rhoe!( Ham, psiks, Rhoe )
        update!( Ham, Rhoe )
        # Calculate energies once again
        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)

        #Etot = Etot_innerscf

        diffE = abs( Etot - Etot_old )
        println("")
        @printf("TRDCM: %5d %18.10f %18.10e\n", iter, Etot, diffE)
        println("")

        if diffE < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("TRDCM is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end

        Etot_old = Etot

        # No need to update potential, it is already updated in inner SCF loop
        for i in 1:Nkspin
            if iter > 1
                P[i] = Y[i][:,set4]*G[i][set4,set1]
            else
                P[i] = Y[i][:,set2]*G[i][set2,set1]
            end
        end

        flush(stdout)
    end  # end of DCM iteration
    
    Ham.electrons.ebands = evals

    if verbose && print_final_ebands
        @printf("\n")
        @printf("----------------------------\n")
        @printf("Final Kohn-Sham eigenvalues:\n")
        @printf("----------------------------\n")
        @printf("\n")
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)
    end

    if verbose && print_final_energies
        @printf("\n")
        @printf("-------------------------\n")
        @printf("Final Kohn-Sham energies:\n")
        @printf("-------------------------\n")
        @printf("\n")
        println(Ham.energies)
    end

    return
end



function _dcm_construct_Y_and_G!(
    Ham, iter::Int64,
    dcm_block_sets::DCMBlockSets,
    psiks, R, T, B, P, Y, G
)

    set1 = dcm_block_sets.set1
    set2 = dcm_block_sets.set2
    set3 = dcm_block_sets.set3
    set5 = dcm_block_sets.set5

    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    for ispin in 1:Nspin, ik in 1:Nkpt

        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
            
        #
        # FIXME: preallocate Hpsi?
        #
        Hpsi = op_H( Ham, psiks[i] )
        #
        # Calculate residual
        #
        psiHpsi = psiks[i]' * Hpsi
        psiHpsi = 0.5*( psiHpsi + psiHpsi' )
        R[i] = Hpsi - psiks[i]*psiHpsi
        Kprec_inplace!( ik, Ham.pw, R[i] )

        #println("R' * R = ")
        #display(R[i]' * R[i]); println()

        #
        # Construct subspace
        #
        Y[i][:,set1] = psiks[i]
        Y[i][:,set2] = R[i]
        #
        if iter > 1
            Y[i][:,set3] = P[i]
        end
            
        #
        # Project kinetic and ionic potential
        #
        if iter > 1
            KY = op_K(Ham, Y[i]) + op_V_Ps_loc(Ham, Y[i]) + op_V_Ps_nloc(Ham, Y[i])
            T[i] = Y[i]'*KY
            bb = Y[i]'*Y[i]
            B[i] = 0.5*( bb + bb' )
        else
            # only set5=1:2*Nstates is active for iter=1
            KY = op_K(Ham, Y[i][:,set5]) + op_V_Ps_loc(Ham, Y[i][:,set5]) +
                 op_V_Ps_nloc(Ham, Y[i][:,set5])
            T[i][set5,set5] = Y[i][:,set5]'*KY
            bb = Y[i][set5,set5]' * Y[i][set5,set5]
            B[i][set5,set5] = 0.5*( bb + bb' )
        end

        # XXX: G is always initialized to identity matrix?
        fill!(G[i], 0.0)
        G[i][set1,set1] = Matrix(1.0I, Nstates, Nstates)
    end

    return
end


function _dcm_project_nonlinear_pot!(
    Ham, iter::Int64,
    dcm_block_sets::DCMBlockSets,
    psiks, σ, evals, A, B, C, D, T, Y, G
)

    set1 = dcm_block_sets.set1
    set2 = dcm_block_sets.set2
    set3 = dcm_block_sets.set3
    set5 = dcm_block_sets.set5

    Nocc = Ham.electrons.Nstates_occ
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    # Nocc should be the same as Nstates

    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        #
        # Project Hartree, XC potential
        #
        V_loc = Ham.potentials.Hartree + Ham.potentials.XC[:,ispin]
        #
        if iter > 1
            yy = Y[i]
        else
            yy = Y[i][:,set5]
        end
        VY = op_V_loc( ik, Ham.pw, V_loc, yy )

        if iter > 1
            A[i] = T[i] + yy'*VY
            A[i] = 0.5*( A[i] + A[i]' )
        else
            aa = T[i][set5,set5] + yy'*VY
            A[i] = 0.5*( aa + aa' )
        end
        #
        if iter > 1
            BG = B[i]*G[i][:,1:Nocc]
            C[i] = BG*BG'
            C[i] = 0.5*( C[i] + C[i]' )
        else
            BG = B[i][set5,set5]*G[i][set5,1:Nocc]
            cc = BG*BG'
            C[i][set5,set5] = 0.5*( cc + cc' )
        end

        if iter > 1
            D[:,i], G[i] =
            eigen( A[i] - σ[i]*C[i], B[i] )
        else
            D[set5,i], G[i][set5,set5] =
            eigen( A[i][set5,set5] - σ[i]*C[i][set5,set5], B[i][set5,set5] )
        end

        evals[:,i] = D[1:Nstates,i]
        #println("evals = ", evals)
        #
        # update wavefunction
        #
        if iter > 1
            psiks[i][:,:] = Y[i]*G[i][:,set1]
            #ortho_sqrt!(psiks[i])  # is this necessary ?
        else
            psiks[i][:,:] = Y[i][:,set5]*G[i][set5,set1]
            #ortho_sqrt!(psiks[i])
        end

        #println("A = ")
        #display(A[i][set5,set5]); println()
        #println("B = ")
        #display(B[i][set5,set5]); println()
        #println("C = ")
        #display(C[i][set5,set5]); println()

        #ortho_check(psiks[i])

        #exit()

    end

    return

end