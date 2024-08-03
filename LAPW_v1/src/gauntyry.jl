function gauntyry(l1, l2, l3, m1, m2, m3)
    c1 = 0.5*sqrt(2) # 0.7071067811865475244d0
    if m2 > 0
        if mod(m2,2)==0
            t1 = c1*( gaunt(l1,l2,l3,m1,m2,m3) + gaunt(l1,l2,l3,m1,-m2,m3) )
        else
            t1 = c1*( gaunt(l1,l2,l3,m1,m2,m3) - gaunt(l1,l2,l3,m1,-m2,m3) )
        end
        res = t1 + im*0.0
    elseif m2 < 0
        if mod(m2,2) == 0
            t1 = c1*( gaunt(l1,l2,l3,m1,m2,m3) - gaunt(l1,l2,l3,m1,-m2,m3) )
        else
            t1 = c1*( gaunt(l1,l2,l3,m1,m2,m3) + gaunt(l1,l2,l3,m1,-m2,m3) )
        end
        res = 0.0 - im*t1
    else
        res = gaunt(l1,l2,l3,m1,m2,m3) + im*0.0
    end
    return res
end
