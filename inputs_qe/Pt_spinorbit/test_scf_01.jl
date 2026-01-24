using PWDFT

atoms = Atoms(1, 1, [0.0; 0.0; 0.0;;], [1], ["Pt"], ["Pt"], [-3.71 0.0 -3.71; 0.0 3.71 3.71; 3.71 3.71 0.0], [18.0], [0.0]);
pspfiles = ["../pseudo_LDA\\Pt_LDA_FR.SG15v1.2.UPF"];
ecutwfc = 15.0;
options = HamiltonianOptions(4.0, 1, 4, [4, 4, 4], [0, 0, 0], false, (0, 0, 0), nothing, "", "VWN", false, -1, 26, true, true, 0.005, [0.0], nothing, nothing, true, true);
pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
for isp in 1:atoms.Nspecies
    pspots[isp] = PsPot_UPF(pspfiles[isp]);
end
Ham = Hamiltonian(atoms, pspots, ecutwfc, options);

psiks = rand_BlochWavefunc(Ham);
electrons_scf_G!(
    Ham,
    psiks = psiks,
    NiterMax = 100,
    betamix = 0.1,
    starting_magn = Ham.options.starting_magn
)
