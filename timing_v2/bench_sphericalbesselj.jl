function bench_sphericalbesselj(Ham)
    pspots = Ham.pspots
    isp = 1
    psp = pspots[isp]
    q = 1.0 # arbitrary ?
    l = 6
    besr = zeros(Float64, psp.kkbeta)
    res = @be begin
        for ir in 1:psp.kkbeta # for PAW this can be kkbeta can be quite extended
            besr[ir] = sphericalbesselj(l, q*psp.r[ir])
        end
    end
    display(res)
    return
end

function bench_qe_sph_bes(Ham)
    pspots = Ham.pspots
    isp = 1
    psp = pspots[isp]
    q = 1.0 # arbitrary ?
    l = 6
    besr = zeros(Float64, psp.kkbeta)
    r = @views psp.r[1:psp.kkbeta]
    res = @be qe_sph_bes!(l, q, r, besr)
    display(res)
    return
end
