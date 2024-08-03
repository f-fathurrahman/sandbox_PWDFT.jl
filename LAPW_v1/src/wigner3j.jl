function wigner3j(j1,j2,j3,m1,m2,m3)
    # check input variables    
    if (j1 < 0) || (j2 < 0) || (j3 < 0) || (abs(m1) > j1) || (abs(m2) > j2)  || (abs(m3) > j3)
        error("Invalid arguments")
    end

    if (j1 == 0) && (j2 == 0) && (j3 == 0)
        return 1.0
    end
    
    if (j1 > 50) || (j2 > 50) || (j3 > 50) then
        error("Angular momenta out of range: $j1 $j1 $j3")
    end

    l1 = j2 - j1 + j3
    l2 = j1 - j2 + j3
    l3 = j1 + j2 - j3
    if (m1+m2+m3 != 0) || (l1 < 0) || (l2 < 0) || (l3 < 0)
        return 0.0
    end
    n1 = j1 - m1
    n2 = j2 + m2
    k1 = max(0, n1-l2, n2-l1)
    k2 = min(l3, n1, n2)
    if mod(k1 - j1 + j2 + m3, 2) != 0
        sgn = -1.0
    else
        sgn = 1.0
    end
    ss = 0.0
    for k in k1:k2
        t1 = sgn * factr(l1,l1-n2+k) * factr(l2,l2-n1+k) * factr(l3,l3-k)
        ss += t1/( factnm(k,1) * factnm(n1-k,1) * factnm(n2-k,1))
        sgn *= -1 # flip sign
    end
    t1 = factr(j1+m1,l1) * factr(j2+m2,l2) * factr(j3+m3,l3)
    t1 = t1*factr(j3-m3,1+j1+j2+j3) * factnm(j1-m1,1) * factnm(j2-m2,1)
    return ss*sqrt(t1)
end
