using Printf

function do_loop(lmax::Int64)
    lm1 = 0
    for l in 0:lmax
        lm2 = lm1 + 2*(l+1)
        for m in -l:-1
            lm1 = lm1+1
            lm2 = lm2-1
            @printf("%3d %3d : %3d %3d\n", l, m, lm1, lm2)
            if mod(m,2) != 0
                #rflm[lm1] = -c1*( imag(zflm[lm1]) + imag(zflm[lm2]) )
                #println("Update lm1 <-  lm1 + lm2 (odd m)")
            else
                #println("Update lm1 <-  lm1 - lm2 (even m)")
                #rflm[lm1] =  c1*( imag(zflm[lm2]) - imag(zflm[lm1]) )
            end
        end
        lm1 = lm1 + 1
        lm2 = lm2 - 1
        #rflm[lm1] = real(zflm[lm1])
        for m in 1:l
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            @printf("%3d %3d : %3d %3d\n", l, m, lm1, lm2)
            if mod(m,2) != 0
                #rflm[lm1] = c1*(real(zflm[lm1]) - real(zflm[lm2]))
                #println("Update lm1 <-  lm1 - lm2 (odd m)")
            else
                #rflm[lm1] = c1*(real(zflm[lm1]) + real(zflm[lm2]))
                #println("Update lm1 <-  lm1 + lm2 (even m)")
            end
        end
    end
    return
end

function main()
    lm = 2
    do_loop(lm)
end

main()