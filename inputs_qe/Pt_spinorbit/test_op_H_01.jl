using PWDFT

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");
#Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv_nospinorb");

psiks = rand_BlochWavefunc(Ham);
Rhoe, _ = atomic_rho_g( Ham, starting_magnetization=[0.0], angle1=[0.0], angle2=[0.0] );
update_from_rhoe!(Ham, psiks, Rhoe);

K_psiks = op_K(Ham, psiks);
Vloc_psiks = op_V_loc(Ham, psiks);
