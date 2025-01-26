# Various functions to update Hamiltonian


# Input: ebands
# Modifies: Focc, E_fermi, mTS
# Also set kT (hardcoded)
function update_from_ebands!(Ham, ebands, kT)

    # NOTE: ebands are assumed to be updated outside this function

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    Nelectrons = Ham.electrons.Nelectrons
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk

    E_fermi, mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT,
        Nkpt, wk
    )

    #println("-----------------------")
    #println("In update_from_ebands!:")
    #println("-----------------------")
    #println()
    #println("E_fermi = ", Ham.electrons.E_fermi)
    #println("mTS = ", Ham.energies.mTS)
    #for ist in 1:Nstates
    #    @printf("%3d %18.10f %18.10f\n", ist, Focc[ist,1], ebands[ist,1])
    #end
    #println("sum Focc: ", sum(Focc))
    #println()
    #println("EXIT update_from_ebands!\n")

    return E_fermi, mTS
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