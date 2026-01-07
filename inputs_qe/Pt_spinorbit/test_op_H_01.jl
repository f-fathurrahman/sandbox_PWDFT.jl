using PWDFT

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");
#Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv_nospinorb");

psiks = rand_BlochWavefunc(Ham);
Rhoe, _ = atomic_rho_g( Ham,
    starting_magn = Ham.options.starting_magn,
    angle1 = Ham.options.angle1,
    angle2 = Ham.options.angle2 
);
update_from_rhoe!(Ham, psiks, Rhoe);

K_psiks = op_K(Ham, psiks);
Vloc_psiks = op_V_loc(Ham, psiks);
Vpsnloc_psiks = op_V_Ps_nloc(Ham, psiks);
