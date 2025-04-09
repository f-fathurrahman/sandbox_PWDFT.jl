# Need to run setup first

function main_scf(filename)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename);

    # This will take into account whether the overlap operator is needed or not
    psiks = rand_BlochWavefunc(Ham)

    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
        Ham.electrons.use_smearing = true
        Ham.electrons.kT = kT
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

    #@infiltrate

    return
end

