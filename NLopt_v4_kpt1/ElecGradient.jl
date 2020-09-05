#
# Electronic gradient: w.r.t psiks and Haux
#
mutable struct ElecGradient
    psiks::BlochWavefunc
    Haux::Vector{Matrix{ComplexF64}}
end

function ElecGradient(Ham)
    psiks = zeros_BlochWavefunc(Ham)
    Nkspin = length(psiks)
    Nstates = Ham.electrons.Nstates
    Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Haux[i] = zeros(ComplexF64,Nstates,Nstates)
    end
    return ElecGradient(psiks, Haux)
end

import Base: -
function -(e1::ElecGradient, e2::ElecGradient)
    return ElecGradient( e1.psiks .- e2.psiks, e1.Haux .- e2.Haux )
end

import Base: length
function length(e::ElecGradient)
    return length(e.psiks)
end