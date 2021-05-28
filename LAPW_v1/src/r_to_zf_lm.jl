function r_to_zf_lm!(lmax, rflm, zflm)
    c1 = 0.7071067811865475244
    lm1 = 0
    for l in 0:lmax
        lm2=lm1+2*(l+1)
        for m in -l:-1
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            if mod(m,2) != 0
                zflm[lm1] = c1*( -rflm[lm2] - im*rflm[lm1] )
            else
                zflm[lm1] = c1*( rflm[lm2] - im*rflm[lm1] )
            end
        end
        lm1 = lm1 + 1
        lm2 = lm2 - 1
        zflm[lm1] = rflm[lm1]
        for m in 1:l
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            if mod(m,2) != 0
                zflm[lm1] = c1*( rflm[lm1] - im*rflm[lm2] )
            else
                zflm[lm1] = c1*( rflm[lm1] + im*rflm[lm2] )
            end
        end
    end
    return
end
