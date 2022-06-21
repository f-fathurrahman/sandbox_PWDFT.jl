using LinearAlgebra
using Printf

using Random
Random.seed!(1234)

N = 20

H = rand(ComplexF64,N,N)
H = 0.5*(H + H')
H = H + N*I # Make diagonally dominant

S = rand(ComplexF64,N,N)
S = 0.5*(S + S')
S = S + N*I # Make diagonally dominant

# res = eigen(H, S)

Nvec = 3
evals = zeros(Float64, Nvec)

# Prepare initial eigenvectors
evc = randn(ComplexF64, N, Nvec)
U = inv(sqrt(evc'*S*evc))
evc = evc*U

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

psi[:,1:Nvec] = evc[:,1:Nvec]

Hpsi[:,1:Nvec] = H*psi[:,1:Nvec]
Spsi[:,1:Nvec] = S*psi[:,1:Nvec]

Hc[1:nbase,1:nbase] = psi[:,1:nbase]'*Hpsi[:,1:nbase]
Sc[1:nbase,1:nbase] = psi[:,1:nbase]'*Spsi[:,1:nbase]

# diagonalize the reduced hamiltonian
ew[1:Nvec], vc[1:Nvec,1:Nvec] = eigen(
    Hermitian(Hc[1:Nvec,1:Nvec]),
    Hermitian(Sc[1:Nvec,1:Nvec])
)
evals[1:Nvec] .= ew[1:Nvec]
for ist in 1:Nvec
    @printf("%4d %18.10f\n", ist, evals[ist])
end

EBANDS_THR = 1e-8
maxter = 40
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
                vc[:,np] .= vc[:,ist]
            end
            # for use in g_psi
            ew[nbase+np] = evals[ist]
        end
    end

    nb1 = nbase + 1
    println("nb1 = ", nb1)
    
    # expand the basis set with new basis vectors ( H - e*S )|psi> ...
    #IF ( uspp ) THEN
    #   CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, spsi, &
    #               kdmx, vc, nvecx, ZERO, psi(:,nb1), kdmx )
    psi[:,nb1:nb1+notcnv-1] = Spsi[:,1:nbase]*vc[1:nbase,1:notcnv]
    for np in 1:notcnv
        psi[:,nbase+np] = -ew[nbase+np]*psi[:,nbase+np]
    end


    #CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, hpsi, &
    #            kdmx, vc, nvecx, ONE, psi(1,nb1), kdmx )
    psi[:,nb1:nb1+notcnv-1] .+= Hpsi[:,1:nbase]*vc[1:nbase,1:notcnv]


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
        psi[:,nbase+n] .= psi[:,nbase+n] / sqrt(ew[n])
    end

    # here compute the hpsi and spsi of the new functions
    Hpsi[:,nb1:nb1+notcnv-1] = H*psi[:,nb1:nb1+notcnv-1]
    Spsi[:,nb1:nb1+notcnv-1] = S*psi[:,nb1:nb1+notcnv-1]
    #Hpsi[:,:] = H*psi
    #Spsi[:,:] = S*psi

    #display(Hpsi); println()
    #display(Spsi); println()

    # update the reduced Hamiltonian
    println("nbase = ", nbase)
    println("notcnv = ", notcnv)
    println("nb1 = ", nb1)
    #CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
    #            kdmx, hpsi(:,nb1), kdmx, ZERO, hc(:,nb1), nvecx )
    Hc[1:nbase+notcnv,nb1:nb1+notcnv-1] = psi[:,1:nbase+notcnv]' * Hpsi[:,nb1:nb1+notcnv-1]

    #CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
    #            kdmx, spsi(:,nb1), kdmx, ZERO, sc(:,nb1), nvecx )
    Sc[1:nbase+notcnv,nb1:nb1+notcnv-1] = psi[:,1:nbase+notcnv]' * Spsi[:,nb1:nb1+notcnv-1]

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

    #CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
    # diagonalize the reduced hamiltonian
    #ew[1:Nvec], vc[1:Nvec,1:Nvec] = eigen(
    #    Hermitian(Hc[1:Nvec,1:Nvec]),
    #    Hermitian(Sc[1:Nvec,1:Nvec])
    #)
    println("Before eigen: nbase = ", nbase)
    ew[1:nbase], vc[1:nbase,1:nbase] = eigen(
        Hermitian(Hc[1:nbase,1:nbase]),
        Hermitian(Sc[1:nbase,1:nbase])
    )
    vc[:,Nvec+1:Nvecx] .= 0.0

    println("ew = ", ew[1:Nvec])
    println("evals  = ", evals[1:Nvec])

    # test for convergence
    for i in 1:Nvec
        is_conv[i] = abs(ew[i] - evals[i]) < EBANDS_THR
    end
    println("is_conv = ", is_conv)
    
    notcnv = sum( .!is_conv )

    #WHERE( btype(1:nvec) == 1 )
    #  conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
    #ELSEWHERE
    #  conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
    #END WHERE


    # Assign new eigenvalues
    evals[1:Nvec] .= ew[1:Nvec]


    # if overall convergence has been achieved, or the dimension of
    # the reduced basis set is becoming too large, or in any case if
    # we are at the last iteration refresh the basis set. i.e. replace
    # the first nvec elements with the current estimate of the
    # eigenvectors;  set the basis dimension to nvec.
    if (notcnv == 0) || ((nbase + notcnv) > Nvecx) || ( dav_iter == maxter )

        println("ENTER THE IF: dav_iter = ", dav_iter)

        #CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, &
        #           psi, kdmx, vc, nvecx, ZERO, evc, kdmx )
        evc[:,1:Nvec] = psi[:,1:nbase]*vc[1:nbase,1:Nvec]

        if notcnv == 0
            println("all roots converged: return")
            println("Should break from loop")
            break
        elseif dav_iter == maxter
            println("Last iteration, some roots not converged: return")
            println("Should break from loop")
            break
        end

        println("Refresh psi, Hpsi, and Spsi")
        println("Should not go here if from last iteration of converged")

        # refresh psi, H*psi and S*psi
        psi[:,1:Nvec] = evc[:,1:Nvec]
        
        #CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, spsi, &
        #            kdmx, vc, nvecx, ZERO, psi(:,nvec+1), kdmx )
        psi[:,Nvec+1:2*Nvec] = Spsi[:,1:nbase]*vc[1:nbase,1:Nvec] 
        Spsi[:,1:Nvec] = psi[:,Nvec+1:2*Nvec]


        #CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, hpsi, &
        #          kdmx, vc, nvecx, ZERO, psi(:,nvec+1), kdmx )
        psi[:,Nvec+1:2*Nvec] = Hpsi[:,1:nbase] * vc[1:nbase,1:Nvec]
        Hpsi[:,1:Nvec] = psi[:,Nvec+1:2*Nvec]

        # refresh the reduced hamiltonian 
        nbase = Nvec
        println("nbase = ", nbase)
        Hc[:,1:nbase] .= 0.0 + im*0.0
        Sc[:,1:nbase] .= 0.0 + im*0.0
        vc[:,1:nbase] .= 0.0 + im*0.0
        for n in 1:nbase
            Hc[n,n] = evals[n] + im*0.0
            Sc[n,n] = 1.0 + im*0.0
            vc[n,n] = 1.0 + im*0.0
        end

    end # if

#=
C := alpha*op( A )*op( B ) + beta*C,

(M x N) = (M x K ) (K x N)

[in]    M   

          M is INTEGER
           On entry,  M specifies the number of rows  of the  matrix
           op( A )  and of the  matrix  C.  M must be at least  zero.

[in]    N   

          N is INTEGER
           On entry,  N  specifies the number  of columns of the matrix
           op( B ) and the number of columns of the matrix C. N must be
           at least zero.

[in]    K   
          K is INTEGER
           On entry,  K  specifies  the number of columns of the matrix
           op( A ) and the number of rows of the matrix op( B ). K must
           be at least  zero.
=#

    # Update iteration
    kter = kter + 1
end


println("Out of iteration")
