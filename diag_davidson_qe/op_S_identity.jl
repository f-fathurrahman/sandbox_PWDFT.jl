# For identity matrix S, simply copy psi to Spsi
function op_S!( Ham::Hamiltonian, psi, Spsi )
    @views Spsi[:] = psi[:]
    return
end