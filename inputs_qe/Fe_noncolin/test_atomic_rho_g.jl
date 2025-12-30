using PWDFT

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");

# FIXME: need to explicitly set starting_magnetization, angle1, and angle2 to zeros
# This is the same as domag=false in pw.x
Rhoe, _ = atomic_rho_g( Ham,
    starting_magn = Ham.options.starting_magn,
    angle1 = Ham.options.angle1,
    angle2 = Ham.options.angle2 
);

psiks = rand_BlochWavefunc(Ham);
update_from_rhoe!(Ham, psiks, Rhoe);
