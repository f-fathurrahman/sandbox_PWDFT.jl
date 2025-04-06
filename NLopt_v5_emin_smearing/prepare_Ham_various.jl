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
        Ham.electrons.use_smearing = true
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
        _, _ = PWDFT._prepare_scf!(Ham, psiks)
        # some heuristics
        is_using_gth_analytic = (eltype(Ham.pspots) == PsPot_GTH)
        is_using_gth_numeric = contains(uppercase(Ham.pspots[1].pspfile), "GTH")
        is_using_hgh_numeric = contains(uppercase(Ham.pspots[1].pspfile), "HGH")
        #
        @assert is_using_gth_analytic || is_using_gth_numeric || is_using_hgh_numeric
        #
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

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP_GTH = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_PSP_GBRV_LDA = joinpath(DIR_PWDFT, "pseudopotentials", "GBRV_LDA")
const DIR_PSP_PAW_JTH_LDA = joinpath(DIR_PWDFT, "pseudopotentials", "PAW_JTH_LDA")


function create_Ham_O2_smearing()
    atoms = Atoms(ext_xyz_file=joinpath(DIR_STRUCTURES, "O2.xyz"))
    pspfiles = [joinpath(DIR_PSP, "O-q6.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4, Nspin=2 )
    Ham.electrons.use_smearing = true
    Ham.electrons.kT = 0.003
    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    #
    Rhoe = guess_rhoe_atomic(Ham, starting_magnetization=[0.1])
    psiks = rand_BlochWavefunc(Ham)
    update!(Ham, psiks, Rhoe)
    return Ham
end

function create_Ham_Pt_fcc_smearing(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    pspfiles = [joinpath(DIR_PSP_GTH, "Pt-q10.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=meshk, extra_states=5 )
    Ham.electrons.use_smearing = true
    Ham.electrons.kT = 0.003
    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    #
    Rhoe = guess_rhoe_atomic(Ham)
    psiks = rand_BlochWavefunc(Ham)
    update!(Ham, psiks, Rhoe)
    return Ham
end


function create_Ham_Pt_fcc_smearing_gbrv(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP_GBRV_LDA, "pt_lda_v1.4.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    #
    options = HamiltonianOptions()
    options.extra_states = 4
    options.dual = ecutrho/ecutwfc
    options.meshk = meshk
    #
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    Ham.electrons.use_smearing = true
    Ham.electrons.kT = 0.003
    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    #
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham)
    _, _ = PWDFT._prepare_scf!(Ham, psiks)
    return Ham
end


function create_Ham_Fe_bcc_smearing_gbrv(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Fe  0.0  0.0  0.0
        """, LatVecs=gen_lattice_bcc(2.866*ANG2BOHR))
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP_GBRV_LDA, "fe_lda_v1.5.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    #
    options = HamiltonianOptions()
    options.extra_states = 4
    options.dual = ecutrho/ecutwfc
    options.meshk = meshk
    options.Nspin = 2
    #
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    Ham.electrons.use_smearing = true
    Ham.electrons.kT = 0.003
    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    #
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham)
    _, _ = PWDFT._prepare_scf!(Ham, psiks)
    return Ham
end


function create_Ham_Fe_bcc_smearing_paw_jth(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Fe  0.0  0.0  0.0
        """, LatVecs=gen_lattice_bcc(2.866*ANG2BOHR))
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP_PAW_JTH_LDA, "Fe.upf"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    #
    options = HamiltonianOptions()
    options.extra_states = 4
    options.dual = ecutrho/ecutwfc
    options.meshk = meshk
    options.Nspin = 2
    #
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    Ham.electrons.use_smearing = true
    Ham.electrons.kT = 0.003
    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    #
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham)
    _, _ = PWDFT._prepare_scf!(Ham, psiks)
    return Ham
end
