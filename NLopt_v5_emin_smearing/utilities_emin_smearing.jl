function get_diag_Haux_from_ebands( Ham )
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
    end
    return Haux
end

function transform_psiks_Haux!(Ham, psiks, Haux)
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    ebands = Ham.electrons.ebands
    Urot = zeros(Float64, Nstates, Nstates)
    for ikspin in 1:Nkspin
        ebands[:,ikspin], Urot[:,:] = eigen(Hermitian(Haux[ikspin]))
        psiks[ikspin] *= Urot
        Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
    end
    return
end

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