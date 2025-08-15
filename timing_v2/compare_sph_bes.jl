function compare_sph_bes(Ham; q=1.0, l=0)

    # the differences are more pronounced for large l value

    pspots = Ham.pspots
    isp = 1
    psp = pspots[isp]
    
    besr1 = zeros(Float64, psp.kkbeta)
    for ir in 1:psp.kkbeta # for PAW this can be kkbeta can be quite extended
        besr1[ir] = sphericalbesselj(l, q*psp.r[ir])
    end
    
    besr2 = zeros(Float64, psp.kkbeta)
    r = @views psp.r[1:psp.kkbeta]
    qe_sph_bes!(l, q, r, besr2)

    println("max abs diff = ", maximum(abs.(besr1-besr2)))
    println("mean abs diff = ", mean(abs.(besr1-besr2)))

    @exfiltrate # to export some variables for plotting

    return
end
