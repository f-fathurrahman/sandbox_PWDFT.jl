function zpotcoul!(
    atoms, atsp_vars,
    mt_vars, pw, 
    zrhoir, # input
    zvclmt, # input output
    zvclir 
)


    println("enter zpotcoul")

    Nspecies = atoms.Nspecies
    rmt = mt_vars.rmt
    gmaxvr = sqrt(2*pw.ecutrho)
    lmaxo = mt_vars.lmaxo

    # Poisson solver pseudocharge density constant
    if Nspecies > 0
        t1 = 0.250*gmaxvr*maximum(rmt)
    else
        t1 = 0.25*gmaxvr*2.0  # FIXME: for which case is this useful?
    end
    println("t1 = ", t1)
    npsd = max(round(Int64,t1),1)
    lnpsd = lmaxo + npsd + 1

    println("npsd = ", npsd)
    println("lnpsd = ", lnpsd)

    return

end