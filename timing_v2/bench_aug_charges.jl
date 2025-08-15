function bench_aug_charges(Ham)
    atoms = Ham.atoms;
    pw = Ham.pw;
    pspots = Ham.pspots;
    #
    pspotNL = Ham.pspotNL;
    lmaxkb = pspotNL.lmaxkb;
    nhm = pspotNL.nhm;
    nh = pspotNL.nh;
    indv = pspotNL.indv;
    nhtolm = pspotNL.nhtolm
    lpl = pspotNL.lpl;
    lpx = pspotNL.lpx;
    ap = pspotNL.ap;
    #
    println("lmaxkb = ", lmaxkb)
    println("nhm = ", nhm)
    println("size ap = ", size(ap))
    println("size indv = ", size(indv))
    println("size nh = ", size(nh))
    println("nh = ", nh)
    #
    res = @be PWDFT._prepare_aug_charges(
        atoms, pw, pspots, lmaxkb, nhm, nh, indv, nhtolm, lpl, lpx, ap
    )
    display(res)
    return
end