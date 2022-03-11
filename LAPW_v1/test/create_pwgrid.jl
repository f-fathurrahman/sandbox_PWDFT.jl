function create_pwgrid(
    atoms, sym_info, mt_vars; rgkmax=7.0, gmaxvr=12.0, kpt_grid=[1,1,1]
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # use average muffin-tin radius (default)
    println("Using average muffin-tin radius to determine gkmax")
    rsum = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        rsum = rsum + mt_vars.rmt[isp]
        println("rmt[isp] = ", mt_vars.rmt[isp])
    end
    rsum = rsum/Natoms
    gkmax = rgkmax/rsum
    println("gkmax = ", gkmax)

    if gmaxvr <= 2*gkmax
        println("gkmax is larger than given/default gmaxvr")
        println("Using gmaxvr = 2*gkmax")
        gmaxvr = 2*gkmax
    end
    println("gmaxvr = ", gmaxvr)

    ecutrho = 0.5*gmaxvr^2
    ecutwfc = 0.5*gkmax^2

    dual = ecutrho/ecutwfc
    println("dual = ", dual)

    pw = PWGrid( ecutwfc, atoms.LatVecs, dual=dual,
        kpoints=KPoints(atoms, kpt_grid, [0,0,0], sym_info.s)
    )
    return pw
end