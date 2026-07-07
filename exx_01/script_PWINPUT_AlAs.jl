atoms_tuple = (Natoms = 2, Nspecies = 2, positions = [2.6550652062682487 0.5310130412536498; 2.6550652062682487 0.5310130412536498; 2.6550652062682487 0.0], atm2species = [1, 2], atsymbs = ["As", "Al"], SpeciesSymbols = ["As", "Al"], LatVecs = [5.310130412536497 5.310130412536497 0.0; 5.310130412536497 0.0 5.310130412536497; 0.0 5.310130412536497 5.310130412536497], Zvals = [15.0, 3.0], masses = [0.0, 0.0]);
atoms = Atoms(atoms_tuple...);
pspfiles = ["./As.upf", "./Al.upf"];
ecutwfc = 20.0;
options_tuple = (dual = 4.0, Nspin_wf = 1, Nspin_dens = 1, meshk = [3, 3, 3], shiftk = [0, 0, 0], time_reversal = true, Ns = (0, 0, 0), kpoints = nothing, kpts_str = nothing, xcfunc = "VWN", use_xc_internal = false, extra_states = nothing, Nstates = nothing, use_symmetry = true, use_smearing = true, smearing_kT = 0.005, starting_magn = nothing, angle1 = nothing, angle2 = nothing, lspinorb = false, noncollinear = false);
options = HamiltonianOptions(options_tuple...);
pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
for isp in 1:atoms.Nspecies
    pspots[isp] = PsPot_UPF(pspfiles[isp]);
end
Ham = Hamiltonian(atoms, pspots, ecutwfc, options);

