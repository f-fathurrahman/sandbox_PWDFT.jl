function print_vec_mat( v::Vector{Matrix{ComplexF64}} )
    Nkspin = length(v)
    for i in 1:Nkspin
        println("Real part of ikspin = ", i)
        display(real(v[i])); println()
        println("Imag part of ikspin = ", i)
        display(imag(v[i])); println()    
    end
end

# 
# Store psiks and Hsub
#
mutable struct ElecVars
    psiks::BlochWavefunc
    Hsub::Array{Matrix{ComplexF64},1}
    Hsub_eigs::Array{Float64,2}
end

function ElecVars( Ham::Hamiltonian )
    return ElecVars( Ham, rand_BlochWavefunc(Ham) )
end

function ElecVars( Ham::Hamiltonian, psiks::BlochWavefunc )
    
    Nkspin = length(psiks)
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin

    Hsub = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    Hsub_eigs = zeros(Float64,Nstates,Nkspin) # the same as electrons.ebands
    
    Rhoe = calc_rhoe( Ham, psiks )
    update!(Ham, Rhoe)

    for ispin in 1:Nspin, ik in 1:Nkpt
        i = ik + (ispin - 1)*Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        #
        Hsub[i] = zeros(ComplexF64,Nstates,Nstates)
        #
        Hsub[i][:] = psiks[i]' * op_H(Ham, psiks[i])
        #
        Hsub_eigs[:,i] = eigvals(Hermitian(Hsub[i]))  # set Haux_eigs to eigenvalues of Hsub
    end

    return ElecVars(psiks, Hsub, Hsub_eigs)
end

function calc_Hsub_eigs!( evars::ElecVars )
    Nkspin = length(evars.Hsub)
    for i in 1:Nkspin
        evars.Hsub_eigs[:,i] = eigvals(Hermitian(evars.Hsub[i]))
    end
end

import Base: show
function show( io::IO, evars::ElecVars )
    Nkspin = length(evars.Hsub)
    Nstates = size(evars.Hsub[1],1)
    for i in 1:Nkspin
        println("Hsub_eigs i = ", i)
        for ist in 1:Nstates
            @printf("%3d %18.10f\n", ist, evars.Hsub_eigs[ist,i])
        end
    end
end
show( evars::ElecVars ) = show( stdout, evars )

# Not yet adapted for Nspin=2
function print_ebands_Hsub_eigs( Ham, evars )
    ebands = Ham.electrons.ebands
    Hsub_eigs = evars.Hsub_eigs
    Nstates, Nkspin = size(ebands)
    println("\nebands and Hsub eigenvalues:")
    for i in 1:Nkspin
        println("ikspin = ", i)
        for ist in 1:Nstates
            @printf("%4d %10.5f %18.10f %18.10f %18.10e\n", ist, Ham.electrons.Focc[ist,i],
                ebands[ist,i], Hsub_eigs[ist,i], abs(ebands[ist,i]-Hsub_eigs[ist,i]))
        end
    end
end

function dot_ElecGradient( v1::ElecGradient, v2::ElecGradient )
    Nkspin = length(v1.psiks)
    ss = 0.0
    for i in 1:Nkspin
        ss = ss + 2.0*real( dot(v1.psiks[i], v2.psiks[i]) )
        ss = ss + real( dot(v1.Haux[i], v2.Haux[i]) ) # no factor of 2
    end
    return ss
end

function dot_ElecGradient_v2( v1::ElecGradient, v2::ElecGradient )
    Nkspin = length(v1.psiks)
    ss = 0.0
    ss_Haux = 0.0
    for i in 1:Nkspin
        ss = ss + 2.0*real( dot(v1.psiks[i], v2.psiks[i]) )
        ss_Haux = ss_Haux + real( dot(v1.Haux[i], v2.Haux[i]) ) # no factor of 2
    end
    return ss, ss_Haux
end