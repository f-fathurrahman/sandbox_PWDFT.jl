using Revise, Infiltrator
using PWDFT

function main_scf()
    #Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");
    #Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_Al");
    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_Al_nospin");

    # This will take into account whether the overlap operator is needed or not
    psiks = rand_BlochWavefunc(Ham)

    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
    end

    if pwinput.nspin == 2
        starting_magnetization = pwinput.starting_magnetization
    else
        starting_magnetization = nothing
    end

    PWDFT.electrons_scf_G!(
        Ham, psiks,
        NiterMax=100,
        use_smearing=use_smearing,
        kT=kT,
        betamix=0.1,
        starting_magnetization=starting_magnetization
    )

    @infiltrate

    return
end

