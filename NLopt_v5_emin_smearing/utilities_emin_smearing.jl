# https://discourse.julialang.org/t/how-to-generate-a-random-unitary-matrix-perfectly-in-julia/34102
function random_unitary_matrix(T, N::Int64)
    #X = (rand(Float64, N, N) + rand(Float64, N, N)*im) / sqrt(2)
    X = rand(Float64, N,N)
    F = qr(X)
    diagR = sign.(real(diag(F.R)))
    diagR[diagR.==0] .= 1
    diagRm = diagm(diagR)
    U = F.Q * diagRm
    return U
end