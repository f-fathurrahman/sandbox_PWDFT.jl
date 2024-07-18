function ylmroty!(β::Float64, lmax, dy)

    ld = size(dy, 1)

    cb = cos(β/2)
    sb = sin(β/2)
    lm1 = 0
    for l in 0:lmax
        # generate rotation operator for m-components of current l
        for m1 in -l:l
            lm1 = lm1 + 1
            lm2 = l^2
            for m2 in -l:l
                lm2 = lm2 + 1
                ss = 0.0
                for k in 0:min(l+m1,l-m2)
                    if ((l+m1-k) >= 0)  &&  ((l-m2-k) >= 0)  &&  ((m2-m1+k) >= 0)
                        j = 2*(l-k) + m1 - m2
                        if j == 0
                            t1 = 1.0
                        else
                            t1 = cb^j
                        end
                        j = 2*k + m2 - m1
                        if j != 0
                            t1 = t1*sb^j
                        end
                        t2 = t1/(factnm(k,1)*factnm(l+m1-k,1)*factnm(l-m2-k,1) *factnm(m2-m1+k,1))
                        if mod(k,2) != 0
                            t2 = -t2
                        end
                        ss += t2
                    end # if
                end # k
                t1 = sqrt(factnm(l+m1,1)*factnm(l-m1,1)*factnm(l+m2,1)*factnm(l-m2,1))
                dy[lm1,lm2] = t1*ss
            end # for m2
        end # for m1
    end
    return
end
