using Printf

function main()

    lmax = 2
    lmaxi = 1
    lmmaxi = (lmaxi+1)^2
    lmaxo = 6
    lmmaxo = (lmaxo+1)^2

    nrci = 70
    nrc = 100
    irco = nrci + 1

    lm = 0
    for l in 0:lmax
        #z2 = im^l * z1
        for m in -l:l
            lm = lm + 1
            @printf("l, m, lm = %3d %3d %3d\n", l, m, lm)
            #z3 = z2*conj(ylmg[lm,ig])
            i = lm
            for irc in 1:nrci
                #zfmt[i] = zfmt[i] + jl[l,irc]*z3
                i = i + lmmaxi
                #@printf("%5d %5d\n", irc, i)
            end
            for irc in irco:nrc
                #zfmt[i] = zfmt[i] + jl[l,irc]*z3
                i = i + lmmaxo
                #@printf("%5d %5d\n", irc, i)
            end
            #println("i = ", i - lm) # minus the offset lm
        end
    end
end

main()