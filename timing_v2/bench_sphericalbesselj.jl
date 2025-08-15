function bench_sphericalbesselj(Ham)
    pspots = Ham.pspots
    isp = 1
    psp = pspots[isp]
    q = 1.0 # arbitrary ?
    l = 1
    besr = zeros(Float64, psp.kkbeta)
    res = @be begin
        for ir in 1:psp.kkbeta # for PAW this can be kkbeta can be quite extended
            besr[ir] = sphericalbesselj(l, q*psp.r[ir])
        end
    end
    display(res)
    return
end
