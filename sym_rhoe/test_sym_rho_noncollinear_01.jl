using PWDFT

atoms = Atoms(1, 1, [0.0; 0.0; 0.0;;], [1], ["Fe"], ["Fe"],
    [2.6085 -2.6085 -2.6085; 2.6085 2.6085 -2.6085; 2.6085 2.6085 2.6085],
    [16.0], [0.0]);
pspfiles = ["../inputs_qe/pseudo_LDA\\Fe_oncv.UPF"];
ecutwfc = 25.0;
options = HamiltonianOptions(4.0, 1, 4, [6, 6, 6],
    [0, 0, 0], false, (0, 0, 0), nothing, "", "VWN", false, -1, 24, true,
    false, [0.5], [90.0], [0.0], false, true);
pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
for isp in 1:atoms.Nspecies
    pspots[isp] = PsPot_UPF(pspfiles[isp]);
end
Ham = Hamiltonian(atoms, pspots, ecutwfc, options);
Ham.rhoe[:,:], _ = atomic_rho_g( Ham,
    starting_magn = Ham.options.starting_magn,
    angle1 = Ham.options.angle1,
    angle2 = Ham.options.angle2 
);

Rhoe = Ham.rhoe;
symmetrize_rhoe_v2!( Ham.pw, Ham.sym_info, Ham.rhoe_symmetrizer, Rhoe )