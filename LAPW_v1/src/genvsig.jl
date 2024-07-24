function genvsig!(pw, vsir, cfunir, vsig)
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    Ng = pw.gvec.Ng
    #
    zfft = zeros(ComplexF64, Npoints)
    # multiply potential by characteristic function in real-space
    zfft[:] .= vsir[:] .* cfunir[:]
    # Fourier transform to G-space
    R_to_G!(pw, zfft)
    # store in global array
    for ig in 1:Ng
        ip = idx_g2r[ig]
        vsig[ig] = zfft[ip]
    end
    # Scale with 1/Npoints
    vsig[:] .*= (1/Npoints)

    return
end

