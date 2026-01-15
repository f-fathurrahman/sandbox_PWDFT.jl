using PWDFT

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");

if pwinput.occupations == "smearing"
    Ham.electrons.use_smearing = true
    Ham.electrons.kT = pwinput.degauss*0.5 # convert from Ry to Ha
end

if pwinput.nspin == 2
    starting_magn = pwinput.starting_magnetization
else
    starting_magn = nothing
end

psiks = rand_BlochWavefunc(Ham);
electrons_scf_G!(
    Ham,
    psiks = psiks,
    NiterMax = 100,
    betamix = 0.1,
    starting_magn = starting_magn
)
