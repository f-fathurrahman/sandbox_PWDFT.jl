mutable struct SubspaceRotations
    prev::Vector{Matrix{ComplexF64}}
    prevC::Vector{Matrix{ComplexF64}}
    prevCinv::Vector{Matrix{ComplexF64}}
end

function SubspaceRotations(Nkspin::Int64, Nstates::Int64)
    prev = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    prevC = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    prevCinv = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    for i in 1:Nkspin
        prev[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        prevC[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        prevCinv[i] = diagm( 0 => ones(ComplexF64,Nstates) )
    end
    return SubspaceRotations(prev, prevC, prevCinv)
end