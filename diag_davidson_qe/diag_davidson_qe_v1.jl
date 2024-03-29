#=
Intended to be used with dense matrix for simulating the actual Hamiltonian
and overlap operator.
The preconditioners are not yet effective for this case.
=#

function diag_davidson_qe!(
    H, S, N, Nvec, evals, evc;
    EBANDS_THR=1e-8, maxter=40
)

    Nvecx = 2*Nvec
    psi = zeros(ComplexF64, N, Nvecx)
    Hpsi = zeros(ComplexF64, N, Nvecx)
    Spsi = zeros(ComplexF64, N, Nvecx)

    Sc = zeros(ComplexF64, Nvecx, Nvecx)
    Hc = zeros(ComplexF64, Nvecx, Nvecx)
    vc = zeros(ComplexF64, Nvecx, Nvecx)

    # eigenvalues of reduced Hamiltonian
    ew = zeros(Float64, Nvecx)

    is_conv = zeros(Bool, Nvec)

    notcnv = Nvec
    nbase  = Nvec

    @views psi[:,1:Nvec] = evc[:,1:Nvec]
    @views Hpsi[:,1:Nvec] = H*psi[:,1:Nvec]
    @views Spsi[:,1:Nvec] = S*psi[:,1:Nvec]

    @views Hc[1:nbase,1:nbase] = psi[:,1:nbase]'*Hpsi[:,1:nbase]
    @views Sc[1:nbase,1:nbase] = psi[:,1:nbase]'*Spsi[:,1:nbase]

    # diagonalize the reduced hamiltonian
    ew[1:Nvec], vc[1:Nvec,1:Nvec] = eigen(
        Hermitian(Hc[1:Nvec,1:Nvec]),
        Hermitian(Sc[1:Nvec,1:Nvec])
    )
    @views evals[1:Nvec] .= ew[1:Nvec]


    dav_iter = 0
    kter = 1
    while (kter <= maxter) || (notcnv == 0)

        dav_iter = kter
        println("kter = ", kter)

        np = 0
        for ist in 1:Nvec
            if !is_conv[ist]
                # this root not yet converged ... 
                np = np + 1
                # reorder eigenvectors so that coefficients for unconverged
                # roots come first. This allows to use quick matrix-matrix 
                # multiplications to set a new basis vector (see below)
                if np != ist
                    @views vc[:,np] .= vc[:,ist]
                end
                # for use in g_psi
                ew[nbase+np] = evals[ist]
            end
        end

        nb1 = nbase + 1
    
        # expand the basis set with new basis vectors ( H - e*S )|psi> ...
        @views psi[:,nb1:nb1+notcnv-1] = Spsi[:,1:nbase]*vc[1:nbase,1:notcnv]
        for np in 1:notcnv
            @views psi[:,nbase+np] = -ew[nbase+np]*psi[:,nbase+np]
        end
        @views psi[:,nb1:nb1+notcnv-1] .+= Hpsi[:,1:nbase]*vc[1:nbase,1:notcnv]

        
        #=
        println("psi before = ")
        display(abs.(psi)); println()

        println("sum psi before = ", sum(psi))
        println("nb1 before precond = ", nb1)
        @views _precond_H_diag_new!(psi[:,nb1:nb1+notcnv-1], H, S, ew[nb1:nb1+notcnv-1])
        #@views _precond_H_diag!(psi[:,nb1:nb1+notcnv-1], H, S, ew[nb1:nb1+notcnv-1])
        #@views _precond_H_inv!(psi[:,nb1:nb1+notcnv-1], H, S, ew[nb1:nb1+notcnv-1])
        println("sum psi after = ", sum(psi))

        println("psi after = ")
        display(abs.(psi)); println()
        =#


        # "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
        # order to improve numerical stability of subspace diagonalization
        # (cdiaghg) ew is used as work array :
        #
        #      ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
        for n in 1:notcnv
            nbn = nbase + n
            @views ew[n] = real(dot(psi[:,nbn], psi[:,nbn]))
        end

        for n in 1:notcnv
            @views psi[:,nbase+n] .= psi[:,nbase+n] / sqrt(ew[n])
        end

        # here compute the hpsi and spsi of the new functions
        @views Hpsi[:,nb1:nb1+notcnv-1] = H*psi[:,nb1:nb1+notcnv-1]
        @views Spsi[:,nb1:nb1+notcnv-1] = S*psi[:,nb1:nb1+notcnv-1]

        # update the reduced Hamiltonian
        @views Hc[1:nbase+notcnv,nb1:nb1+notcnv-1] = psi[:,1:nbase+notcnv]' * Hpsi[:,nb1:nb1+notcnv-1]
        @views Sc[1:nbase+notcnv,nb1:nb1+notcnv-1] = psi[:,1:nbase+notcnv]' * Spsi[:,nb1:nb1+notcnv-1]

        nbase = nbase + notcnv
        for n in 1:nbase
            # the diagonal of hc and sc must be strictly real 
            Hc[n,n] = real(Hc[n,n]) + im*0.0
            Sc[n,n] = real(Sc[n,n]) + im*0.0
            #
            for m in n+1:nbase
                Hc[m,n] = conj(Hc[n,m])
                Sc[m,n] = conj(Sc[n,m])
            end
        end

        println("Hc = ")
        display(abs.(Hc[1:nbase,1:nbase])); println()
        println("Sc = ")
        display(abs.(Sc[1:nbase,1:nbase])); println()

        # diagonalize the reduced hamiltonian
        ew[1:nbase], vc[1:nbase,1:nbase] = eigen(
            Hermitian(Hc[1:nbase,1:nbase]),
            Hermitian(Sc[1:nbase,1:nbase])
        )
        @views vc[:,Nvec+1:Nvecx] .= 0.0

        # test for convergence
        # FIXME: Use different EBANDS for occupied and unoccupied bands
        println("Convergence")
        for i in 1:Nvec
            Δ = abs(ew[i] - evals[i])
            @printf("%3d %18.10f %18.10f %18.10e\n", i, ew[i], evals[i], Δ)
            is_conv[i] = abs(ew[i] - evals[i]) < EBANDS_THR
        end
        println("is_conv = ", is_conv)
    

        notcnv = sum( .!is_conv )

        # Assign new eigenvalues
        @views evals[1:Nvec] .= ew[1:Nvec]

        # if overall convergence has been achieved, or the dimension of
        # the reduced basis set is becoming too large, or in any case if
        # we are at the last iteration refresh the basis set. i.e. replace
        # the first nvec elements with the current estimate of the
        # eigenvectors;  set the basis dimension to nvec.
        if (notcnv == 0) || ((nbase + notcnv) > Nvecx) || ( dav_iter == maxter )

            @views evc[:,1:Nvec] = psi[:,1:nbase]*vc[1:nbase,1:Nvec]

            if notcnv == 0
                println("all roots converged in kter = ", dav_iter)
                break
            elseif dav_iter == maxter
                println("Last iteration, some roots not converged: return")
                break
            end

            # refresh psi, H*psi and S*psi
            @views psi[:,1:Nvec] = evc[:,1:Nvec]
            @views psi[:,Nvec+1:2*Nvec] = Spsi[:,1:nbase]*vc[1:nbase,1:Nvec] 
            @views Spsi[:,1:Nvec] = psi[:,Nvec+1:2*Nvec]
            @views psi[:,Nvec+1:2*Nvec] = Hpsi[:,1:nbase] * vc[1:nbase,1:Nvec]
            @views Hpsi[:,1:Nvec] = psi[:,Nvec+1:2*Nvec]

            # refresh the reduced hamiltonian 
            nbase = Nvec
            @views Hc[:,1:nbase] .= 0.0 + im*0.0
            @views Sc[:,1:nbase] .= 0.0 + im*0.0
            @views vc[:,1:nbase] .= 0.0 + im*0.0
            for n in 1:nbase
                Hc[n,n] = evals[n] + im*0.0
                Sc[n,n] = 1.0 + im*0.0
                vc[n,n] = 1.0 + im*0.0
            end

        end # if

        # Update iteration counter
        kter = kter + 1
    end

    return

end



function _precond_H_diag!(psi, H, S, ew)
    SMALL = 1e-4
    Nbasis = size(psi,1)
    println("ew in _precond_H_diag = ", ew)
    notcnv = size(ew,1)

    #scal = 1.0
    #for k in 1:notcnv, i in 1:Nbasis
    #    x = ( H[i,i] - ew[k]*S[i,i] )*scal
    #    denm = 0.5*( 1.0 + x + sqrt(1.0 + (x-1)*(x-1)) )/scal
    #    psi[i,k] = psi[i,k]/denm
    #end

    for ist in 1:notcnv, i in 1:Nbasis
        # denm = g2 + v(g=0) - e(k)
        denm = real(H[i,i] - ew[ist]*S[i,i])
        if abs(denm) < SMALL
            denm = sign(denm)*SMALL
        end
        psi[i,ist] = psi[i,ist]/denm
    end

    return
end


function _precond_H_diag_new!(psi, H, S, ew)
    SMALL = 1e-4
    Nbasis = size(psi,1)
    println("ew in _precond_H_diag = ", ew)
    notcnv = size(ew,1)
    scal = 1.0
    for k in 1:notcnv, i in 1:Nbasis
        x = ( H[i,i] - ew[k]*S[i,i] )*scal
        denm = 0.5*( 1.0 + x + sqrt(1.0 + (x-1)*(x-1)) )/scal
        psi[i,k] = psi[i,k]/denm
    end
    return
end


# This is expensive
function _precond_H_inv!(psi, H, S, ew)
    SMALL = 1e-4
    Nbasis = size(psi,1)
    notcnv = size(ew,1)
    println("notcnv = ", notcnv)
    println("ew in _precond_H_inv: ", ew)
    Hinv = zeros(ComplexF64,Nbasis,Nbasis)
    for ist in 1:notcnv
        @views Hinv[:,:] = Hermitian(inv(H - ew[ist]*S)) # need Hermitian?
        @views psi[:,ist] = Hinv*psi[:,ist]
    end
    return
end