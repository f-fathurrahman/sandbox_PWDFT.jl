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

Nstates = 3
evals = zeros(Float64, Nstates)

# Prepare initial eigenvectors
evc = randn(ComplexF64, N, Nstates)
U = inv(sqrt(evc'*S*evc))
evc = evc*U

Nvecx = 2*Nstates
psi = zeros(ComplexF64, N, Nvecx)
Hpsi = zeros(ComplexF64, N, Nvecx)
Spsi = zeros(ComplexF64, N, Nvecx)

Sc = zeros(ComplexF64, Nvecx, Nvecx)
Hc = zeros(ComplexF64, Nvecx, Nvecx)
vc = zeros(ComplexF64, Nvecx, Nvecx)

# eigenvalues of reduced Hamiltonian
ew = zeros(Float64, Nvecx)

is_conv = zeros(Bool, Nstates)

notcnv = Nstates
nbase  = Nstates


psi[:,1:Nstates] = evc[:,1:Nstates]

Hpsi[:,1:Nstates] = H*psi[:,1:Nstates]
Spsi[:,1:Nstates] = S*psi[:,1:Nstates]

Hc = psi'*Hpsi
Sc = psi'*Spsi

# diagonalize the reduced hamiltonian
ew[1:Nstates], vc[1:Nstates,1:Nstates] = eigen(
    Hermitian(Hc[1:Nstates,1:Nstates]),
    Hermitian(Sc[1:Nstates,1:Nstates])
)
evals[1:Nstates] .= ew[1:Nstates]
for ist in 1:Nstates
    @printf("%4d %18.10f\n", ist, evals[ist])
end


maxter = 1
dav_iter = 0
kter = 1
while (kter <= maxter) || (notcnv == 0)

    dav_iter = kter

    println("kter = ", kter)
    np = 0
    for ist in 1:Nstates
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
    
    psi[:,nb1:nb1+notcnv-1] = Spsi[:,1:nbase]*vc[1:nbase,1:nbase]

    for np in 1:notcnv
        psi[:,nbase+np] .= -ew[nbase+np]*psi[:,nbase+np]
    end


    #CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, hpsi, &
    #            kdmx, vc, nvecx, ONE, psi(1,nb1), kdmx )

    psi[:,nb1:nb1+notcnv-1] .+= Hpsi[:,1:nbase]*vc[1:nbase,1:nbase]
    #psi[:,:] += Hpsi*vc
    println(psi[1,nb1])


    # "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
    # order to improve numerical stability of subspace diagonalization
    # (cdiaghg) ew is used as work array :
    #
    #      ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
    for n in 1:notcnv
        nbn = nbase + n
        #ew[n] = ddot( 2*npw, psi(:,nbn), 1, psi(:,nbn), 1 )
        @views ew[n] = real(dot(psi[:,nbn], psi[:,nbn]))
    end

    for n in 1:notcnv
        psi[:,nbase+n] .= psi[:,nbase+n] / sqrt(ew[n])
    end
    println(psi[1,nb1])


    # here compute the hpsi and spsi of the new functions
    #Hpsi[:,nb1:nb1+notcnv-1] = H*psi[:,nb1:nb1+notcnv-1]
    #Spsi[:,nb1:nb1+notcnv-1] = S*psi[:,nb1:nb1+notcnv-1]
    Hpsi[:,:] = H*psi
    Spsi[:,:] = S*psi

    display(Hpsi); println()
    display(Spsi); println()

    # update the reduced Hamiltonian
    #CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
    #            kdmx, hpsi(:,nb1), kdmx, ZERO, hc(:,nb1), nvecx )
    #Hc[1:notcnv,nb1:nb1+notcnv-1] = psi[:,1:notcnv]' * Hpsi[:,1:notcnv]
    Hc = psi' * Hpsi

    #CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
    #            kdmx, spsi(:,nb1), kdmx, ZERO, sc(:,nb1), nvecx )
    #Sc[1:notcnv,nb1:nb1+notcnv-1] = psi[:,1:notcnv]' * Spsi[:,1:notcnv]
    Sc = psi' * Spsi

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

    #display(Hc); println()
    #display(Sc); println()

    println("Hc[5,5] = ", Hc[5,5])
    println("Hc[6,6] = ", Hc[6,6])

    println("Sc[5,5] = ", Sc[5,5])
    println("Sc[6,6] = ", Sc[6,6])


    #CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
    # diagonalize the reduced hamiltonian
    println("nbase = ", nbase)
    println("Nstates = ", Nstates)
    ew[:], vc[:,:] = eigen(
        Hermitian(Hc),
        Hermitian(Sc)
    )
    println("ew = ", ew)


#=
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


