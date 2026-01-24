using Revise, PWDFT

atoms = Atoms(1, 1, [0.0; 0.0; 0.0;;], [1], ["Fe"], ["Fe"],
    [2.6085 -2.6085 -2.6085; 2.6085 2.6085 -2.6085; 2.6085 2.6085 2.6085], [16.0], [0.0]);
pspfiles = ["../pseudo_LDA\\Fe_oncv.UPF"];
ecutwfc = 25.0;
options = HamiltonianOptions(4.0, 1, 4, [6, 6, 6], [0, 0, 0], false, (0, 0, 0),
    nothing, "", "VWN", false, -1, 24, true, true, 0.005, [0.5], [90.0], [0.0], false, true);
pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
for isp in 1:atoms.Nspecies
    pspots[isp] = PsPot_UPF(pspfiles[isp]);
end
Ham = Hamiltonian(atoms, pspots, ecutwfc, options);

psiks = rand_BlochWavefunc(Ham);
electrons_scf_G!(
    Ham,
    psiks = psiks,
    NiterMax = 50,
    betamix = 0.2,
    starting_magn = Ham.options.starting_magn,
    print_final_ebands = false
)