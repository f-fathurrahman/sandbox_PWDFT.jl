using Revise, PWDFT

filename = "PWINPUT_AlAs";
Ham, pwinput = init_Ham_from_pwinput(filename=filename);
Ham.xc_calc = NoneXCCalculator();

psiks = rand_BlochWavefunc(Ham);
electrons_scf_G!(
    Ham,
    psiks = psiks,
    NiterMax = 100,
    betamix = 0.1,
    starting_magn = Ham.options.starting_magn
)
