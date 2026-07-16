using Printf
using LinearAlgebra: norm, inv, dot
using FFTW
using Serialization: serialize, deserialize

function debug_diag_exx()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    
    Ham.exx = EXXVariables(Ham, pwinput)
    #psiks = deserialize("psiks_nox_noc.jldat")
    psiks = deserialize("psiks_latest.jldat")
    Rhoe = deserialize("Rhoe_latest.jldat")
    #Rhoe = calc_rhoe(Ham, psiks)
    #psiks = deserialize("psiks_nox_noc_v01.jldat")
    Ham.exx.is_active = true
    Ham.xc_calc = NoneXCCalculator() # set XC components to zero
    
    E_fock = Inf
    E_fock_old = E_fock
    dE_fock = Inf
    E_total = sum(Ham.energies)
    E_total_old = E_total
    for iterOuter in 1:10
        #
        set_exx_buffer!(Ham, psiks)
        electrons_scf_G!(
            Ham,
            psiks = psiks,
            Rhoe = Rhoe,
            NiterMax = 100,
            betamix = 0.1,
            starting_magn = Ham.options.starting_magn
        )
        E_fock = calc_E_fock(Ham, psiks)
        E_total = sum(Ham.energies) + E_fock
        #
        dE = abs(E_total - E_total_old)
        dE_fock = abs(E_fock - E_fock_old)
        println()
        println("**** Outer iteration: $iterOuter E_fock = $(E_fock) dE_fock = $(dE_fock)")
        println("                                 E_total = $(E_total) dE = $dE")
        println()
        #
        E_fock_old = E_fock
        E_total_old = E_total
    end

    serialize("psiks_latest.jldat", psiks)
    serialize("Rhoe_latest.jldat", Rhoe)

    @infiltrate
end
