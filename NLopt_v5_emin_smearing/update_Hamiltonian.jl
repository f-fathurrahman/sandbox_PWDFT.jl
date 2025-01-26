# Various functions to update Hamiltonian


# Input: ebands
# Modifies: Focc, E_fermi, mTS
# Also set kT (hardcoded)
function update_from_ebands!(Ham, ebands)

    # NOTE: ebands are assumed to be updated outside this function

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    Nelectrons = Ham.electrons.Nelectrons
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    kT = Ham.electrons.kT
    @assert kT > 1e-10

    E_fermi, mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT,
        Nkpt, wk
    )
    # Set some output
    Ham.electrons.E_fermi = E_fermi
    Ham.energies.mTS = mTS

    return
end


# Input: psiks
# Modifies: Ham.rhoe, potentials
function update_from_wavefunc!(Ham, psiks)    
    # Compute electron density from psiks
    # Use Ham.rhoe
    calc_rhoe!(Ham, psiks, Ham.rhoe)
    # Update the potentials
    update_from_rhoe!(Ham, psiks, Ham.rhoe)
    # XXX: update_from_rhoe! will not overwrite update Ham.rhoe
    return
end