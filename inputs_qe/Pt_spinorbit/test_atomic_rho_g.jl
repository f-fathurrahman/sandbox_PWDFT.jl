using PWDFT

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");

# FIXME: need to explicitly set starting_magnetization, angle1, and angle2 to zeros
# This is the same as domag=false in pw.x
Rhoe, _ = atomic_rho_g( Ham, starting_magnetization=[0.0], angle1=[0.0], angle2=[0.0] );

psiks = rand_BlochWavefunc(Ham)

