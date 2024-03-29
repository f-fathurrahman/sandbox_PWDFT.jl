# !INPUT/OUTPUT PARAMETERS:
#   lmax : maximum angular momentum (in,integer)
#   zflm : coefficients of complex spherical harmonic expansion
#          (in,complex((lmax+1)**2)))
#   rflm : coefficients of real spherical harmonic expansion
#          (out,real((lmax+1)**2)))
function z_to_rf_lm!(lmax::Int64, zflm, rflm)

    c1 = 0.7071067811865475244  # 1/sqrt(2)
    lm1 = 0
    for l in 0:lmax
        lm2 = lm1 + 2*(l+1)
        for m in -l:-1
            lm1 = lm1+1
            lm2 = lm2-1
            if mod(m,2) != 0
                rflm[lm1] = -c1*( imag(zflm[lm1]) + imag(zflm[lm2]) )
            else
                rflm[lm1] =  c1*( imag(zflm[lm2]) - imag(zflm[lm1]) )
            end
        end
        lm1 = lm1 + 1
        lm2 = lm2 - 1
        rflm[lm1] = real(zflm[lm1])
        for m in 1:l
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            if mod(m,2) != 0
                rflm[lm1] = c1*(real(zflm[lm1]) - real(zflm[lm2]))
            else
                rflm[lm1] = c1*(real(zflm[lm1]) + real(zflm[lm2]))
            end
        end
    end
    return
end
