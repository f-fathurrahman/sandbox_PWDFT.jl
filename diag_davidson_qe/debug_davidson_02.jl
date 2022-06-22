using LinearAlgebra
using Printf

using Random

include("diag_davidson_qe_v1.jl")

function main()

    Random.seed!(1234)

    N = 20
    #
    H = rand(ComplexF64,N,N)
    H = 0.5*(H + H')
    H = H + N*I # Make diagonally dominant
    #
    S = rand(ComplexF64,N,N)
    S = 0.5*(S + S')
    S = S + N*I # Make diagonally dominant

    Nvec = 3
    evals = zeros(Float64, Nvec)

    # Prepare initial eigenvectors
    evc = randn(ComplexF64, N, Nvec)
    U = inv(sqrt(evc'*S*evc))
    evc = evc*U

    diag_davidson_qe!(H, S, N, Nvec, evals, evc)

    res = eigen(H, S)
    evals_ref = res.values
    for i in 1:Nvec
        Δ = abs(evals[i] - evals_ref[i])
        @printf("%3d %18.10f %18.10f %18.10e\n", i, evals[i], evals_ref[i], Δ)
    end

end # main

main()
#@time main()
