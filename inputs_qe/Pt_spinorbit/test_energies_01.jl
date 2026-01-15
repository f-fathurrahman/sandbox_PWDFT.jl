using PWDFT

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");
#Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv_nospinorb");

psiks = rand_BlochWavefunc(Ham);
Rhoe = Ham.rhoe; # alias
Rhoe[:,:], _ = atomic_rho_g( Ham,
    starting_magn = Ham.options.starting_magn,
    angle1 = Ham.options.angle1,
    angle2 = Ham.options.angle2 
);
update_from_rhoe!(Ham, psiks, Rhoe);

E_kin = calc_E_kin( Ham, psiks )
E_Ps_loc, E_Hartree, E_xc = calc_E_local( Ham, psiks )
E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )
