#includet("atomic_wfc_01.jl")

function create_Ham_structure_01()
    atoms_tuple = (
        Natoms = 2,
        Nspecies = 1,
        positions = [0.0 2.7079775377810678; 0.0 2.7079775377810678; 0.0 2.7079775377810678],
        atm2species = [1, 1],
        atsymbs = ["Fe", "Fe"],
        SpeciesSymbols = ["Fe"],
        LatVecs = [5.4159550755621355 0.0 0.0; 0.0 5.4159550755621355 0.0; 0.0 0.0 5.4159550755621355],
        Zvals = [16.0],
        masses = [0.0]
    );
    atoms = Atoms(atoms_tuple...);
    
    pspfiles = [ joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "Fe.upf") ];
    ecutwfc = 20.0;
    options_tuple = (
        dual = 4.0, Nspin_wf = 2, Nspin_dens = 2, meshk = [3, 3, 3],
        shiftk = [0, 0, 0],
        time_reversal = true, Ns = (0, 0, 0), kpoints = nothing, kpts_str = nothing,
        xcfunc = "VWN",
        use_xc_internal = false, extra_states = nothing,
        Nstates = 20, use_symmetry = true, use_smearing = true,
        smearing_kT = 0.001, starting_magn = [0.4], angle1 = nothing,
        angle2 = nothing, lspinorb = false, noncollinear = false
    );
    options = HamiltonianOptions(options_tuple...);
    pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
    for isp in 1:atoms.Nspecies
        pspots[isp] = PsPot_UPF(pspfiles[isp]);
    end
    Ham = Hamiltonian(atoms, pspots, ecutwfc, options);
    return Ham
end

# Cu fcc
function create_Ham_structure_02()
    atoms = Atoms( xyz_string_frac=
        """
        1

        Cu  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(3.61496*ANG2BOHR) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "Cu.upf")]
    ecutwfc = 30.0

    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[10,10,10],
        extra_states=4, use_smearing=true, smearing_kT=0.01 )
end

# Au fcc, 4 atoms
function create_Ham_structure_03()
    atoms_tuple = (
        Natoms = 4,
        Nspecies = 1,
        positions = [0.0 0.0 3.9439529101367605 3.9439529101367605;
                     0.0 3.9439529101367605 0.0 3.9439529101367605;
                     0.0 3.9439529101367605 3.9439529101367605 0.0],
        atm2species = [1, 1, 1, 1],
        atsymbs = ["Au", "Au", "Au", "Au"],
        SpeciesSymbols = ["Au"],
        LatVecs = [7.887905820273521 0.0 0.0;
                   0.0 7.887905820273521 0.0;
                   0.0 0.0 7.887905820273521],
        Zvals = [11.0],
        masses = [0.0]
    );
    atoms = Atoms(atoms_tuple...);
    pspfiles = [ joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "Au.upf") ];
    ecutwfc = 20.0;
    options_tuple = (
        dual = 4.0, Nspin_wf = 1, Nspin_dens = 1,
        meshk = [4, 4, 4], shiftk = [0, 0, 0], time_reversal = true,
        Ns = (0, 0, 0),
        kpoints = nothing, kpts_str = nothing,
        xcfunc = "VWN", use_xc_internal = false,
        extra_states = nothing, Nstates = nothing,
        use_symmetry = true, use_smearing = true, smearing_kT = 0.005,
        starting_magn = nothing, angle1 = nothing, angle2 = nothing,
        lspinorb = false, noncollinear = false
    );
    options = HamiltonianOptions(options_tuple...);
    pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
    for isp in 1:atoms.Nspecies
        pspots[isp] = PsPot_UPF(pspfiles[isp]);
    end
    return Hamiltonian(atoms, pspots, ecutwfc, options);
end

# Al fcc, 4 atoms, ONCV pspot
function create_Ham_structure_04()
    atoms_tuple = (
        Natoms = 4,
        Nspecies = 1,
        positions = [0.0 0.0 3.817445194667986 3.817445194667986;
                     0.0 3.817445194667986 0.0 3.817445194667986;
                     0.0 3.817445194667986 3.817445194667986 0.0],
        atm2species = [1, 1, 1, 1],
        atsymbs = ["Al", "Al", "Al", "Al"],
        SpeciesSymbols = ["Al"],
        LatVecs = [7.634890389335972 0.0 0.0;
                   0.0 7.634890389335972 0.0;
                   0.0 0.0 7.634890389335972],
        Zvals = [3.0],
        masses = [0.0]
    );
    atoms = Atoms(atoms_tuple...);
    pspfiles = [ joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "Al.upf") ];
    ecutwfc = 15.0;
    options_tuple = (
        dual = 4.0, Nspin_wf = 1, Nspin_dens = 1,
        meshk = [4, 4, 4], shiftk = [0, 0, 0], time_reversal = true,
        Ns = (0, 0, 0), kpoints = nothing, kpts_str = nothing,
        xcfunc = "VWN", use_xc_internal = false, extra_states = nothing,
        Nstates = nothing, use_symmetry = true,
        use_smearing = true, smearing_kT = 0.005,
        starting_magn = nothing, angle1 = nothing, angle2 = nothing,
        lspinorb = false, noncollinear = false
    );
    options = HamiltonianOptions(options_tuple...);
    pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
    for isp in 1:atoms.Nspecies
        pspots[isp] = PsPot_UPF(pspfiles[isp]);
    end
    return Hamiltonian(atoms, pspots, ecutwfc, options);
end


#=

# Calculate atomic_wfc
Nstates = Ham.electrons.Nstates;
Nspin = Ham.electrons.Nspin_wf;
Nkpt = Ham.pw.gvecw.kpoints.Nkpt;
Natomwfc = calc_Natomwfc(Ham.atoms, Ham.pspots);
println("Natomwfc = ", Natomwfc);
println("Nstates = ", Nstates);
psiks = zeros_BlochWavefunc(Ham);
for ispin in 1:Nspin, ik in 1:Nkpt
    ikspin = ik + (ispin-1)*Nkpt
    atomic_wfc!(ik, Ham.atoms, Ham.pspots, Ham.pw, psiks[ikspin]);
end

# XXX Explicitly diagonalize?
for ispin in 1:Nspin, ik in 1:Nkpt
    ikspin = ik + (ispin-1)*Nkpt
    ortho_sqrt!(Ham, psiks[ikspin])
end

Nkspin = Nkpt*Nspin;
Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin);
for ikspin in 1:Nkspin
    Haux[ikspin] = randn(ComplexF64, Nstates, Nstates);
    Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' );
end

# Cu-fcc
# iterCG: 4 E_new = -189.35662141833518 ΔE = 6.134200702945236e-9



Ham.rhoe[:,:], _ = atomic_rho_g(
    Ham,
    starting_magn = Ham.options.starting_magn,
    angle1 = Ham.options.angle1,
    angle2 = Ham.options.angle2
)
update_from_rhoe!(Ham, psiks, Ham.rhoe)

ik = 1; ispin = 1;
for ispin in 1:Nspin, ik in 1:Nkpt
    Ham.ispin = ispin
    Ham.ik = ik
    ikspin = ik + (ispin-1)*Nkpt
    psi = psiks[ikspin];
    Hpsi = op_H(Ham, psiks[ikspin]);
    Hsub = psi' * Hpsi;
    λ, U = eigen(Hsub);
end
=#