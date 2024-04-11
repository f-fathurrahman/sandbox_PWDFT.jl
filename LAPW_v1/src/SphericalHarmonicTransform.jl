struct SphericalHarmonicTransform
    rbshti::Matrix{Float64}
    rfshti::Matrix{Float64}
    zbshti::Matrix{ComplexF64}
    zfshti::Matrix{ComplexF64}
    rbshto::Matrix{Float64}
    rfshto::Matrix{Float64}
    zbshto::Matrix{ComplexF64}
    zfshto::Matrix{ComplexF64}
    # fixed to false for the "normal" calculation
    # will be used in MAE (magnetic anisotropy energy) calculation
    trotsht::Bool
end

function SphericalHarmonicTransform(lmaxi, lmaxo)
    rbshti, zbshti, rfshti, zfshti = _gen_sht_matrices(lmaxi)
    rbshto, zbshto, rfshto, zfshto = _gen_sht_matrices(lmaxo)
    trotsht = false
    return SphericalHarmonicTransform(
        rbshti, rfshti, zbshti, zfshti,
        rbshto, rfshto, zbshto, zfshto,
        trotsht
    )
end


function _gen_sht_matrices(lmax)
    #
    lmmax = (lmax+1)^2
    tp = zeros(2,lmmax)
    vtp = zeros(3,lmmax)
    # generate spherical covering set
    sphcover!(lmmax, tp)
    # convert (theta, phi) angles to vectors
    sctovec!( lmmax, tp, vtp )

    rlm = zeros(Float64,lmmax)
    ylm = zeros(ComplexF64,lmmax)
    #
    # Backward SHT matrices (complex and real)
    #
    rbsht = zeros(Float64,lmmax,lmmax)
    zbsht = zeros(ComplexF64,lmmax,lmmax)
    for lm in 1:lmmax
        @views genrlmv!(lmax, vtp[:,lm], rlm)
        rbsht[lm,:] .= rlm[:]
        #
        @views genylmv!(lmax, vtp[:,lm], ylm)
        zbsht[lm,:] .= ylm[:]
    end

    rfsht = inv(rbsht)
    zfsht = inv(zbsht)

    return rbsht, zbsht, rfsht, zfsht
end