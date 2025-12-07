function init_rot_ylm( ; lmaxx=3 )
    sqrt2 = sqrt(2.0)
    lqmax = 2*lmaxx + 1
    rot_ylm = zeros(ComplexF64, lqmax, lqmax)
    l = lmaxx
    rot_ylm[l+1,1] = 1.0
    for n1 in range(2, stop=2*l+1, step=2)
        m = floor(Int64, n1/2)
        n = l + 1 - m
        rot_ylm[n,n1] = complex( (-1)^m/sqrt2, 0.0 )
        rot_ylm[n,n1+1] = complex( 0.0, -(-1)^m/sqrt2 )
        n = l + 1 + m
        rot_ylm[n,n1] = complex(1.0/sqrt2, 0.0)
        rot_ylm[n,n1+1] = complex(0.0, 1.0/sqrt2)
    end
    return rot_ylm
end