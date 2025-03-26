function prepare_Ham_from_pwinput(filename::String)

    Ham, pwinput = init_Ham_from_pwinput(filename=filename);

    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms);
    
    # We need to set some parameters manually:
    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
        Ham.electrons.kT = kT
    end

    if pwinput.nspin == 2
        starting_magnetization = pwinput.starting_magnetization
    else
        starting_magnetization = nothing
    end

    Nspin = Ham.electrons.Nspin
    # Initialize electronic variables: `psiks` and `Haux`:
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham)
    # XXX: this is not really needed, should be able to pass nothing to

    # Initial Rhoe, taking into account starting magnetization
    #
    # XXX: This is needed for GTH pspots because rho_atom is not available
    if Nspin == 2
        _, _ = update_from_rhoe!( Ham, psiks, Ham.rhoe )
        Rhoe_tot = Ham.rhoe[:,1] + Ham.rhoe[:,2]
        starting_magnetization = 2.6
        magn = starting_magnetization .* ones(size(Rhoe_tot)) / Ham.pw.CellVolume
        Ham.rhoe[:,1] .= 0.5*(Rhoe_tot + magn)
        Ham.rhoe[:,2] .= 0.5*(Rhoe_tot - magn)
        # Update again the Hamiltonian
        _, _ = update_from_rhoe!( Ham, psiks, Ham.rhoe )
    else
        _, _ = PWDFT._prepare_scf!(Ham, psiks)
    end

    return Ham
end