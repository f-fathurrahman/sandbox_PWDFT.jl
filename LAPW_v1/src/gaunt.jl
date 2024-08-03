function gaunt(l1, l2, l3, m1, m2, m3)
    c1 = 0.5/sqrt(pi)
    # XXX more details
    if (l1 < 0) || (l2 < 0) || (l3 < 0) || (abs(m1) > l1) || (abs(m2) > l2) || (abs(m3) > l3)
        error("Error(gaunt): non-physical arguments")
    end
    #
    if (l1 > 50) || (l2 > 50) || (l3 > 50)
        error("Error(gaunt): angular momenta out of range: $l1 $l2 $l3")
    end
    #
    if m1-m2-m3 != 0
        return 0.0
    end
    j1 = l2 - l1 + l3
    j2 = l1 - l2 + l3
    j3 = l1 + l2 - l3
    #
    if (j1 < 0) || (j2 < 0) || (j3 < 0)
        return 0.0
    end
    j = l1 + l2 + l3
    if mod(j,2) != 0
        return 0.0
    end
    jh = j/2
    t1 = sqrt( (2*l1+1)*(2*l2+1)*(2*l3+1)*factr(j1,j+1)*factnm(j2,1)*factnm(j3,1) )
    t1 = t1*factr(jh,jh-l1) / ( factnm(jh-l2,1)*factnm(jh-l3,1) )
    res = t1*c1*wigner3j(l1, l2, l3, -m1, m2, m3)
    if mod(m1+jh,2) != 0
        res *= -1
    end
    return res
end
