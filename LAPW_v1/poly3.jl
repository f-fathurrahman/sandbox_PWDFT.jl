function poly3(xa, ya, x)
    # arguments
    # real(8) xa(3),ya(3),x
    # local variables

    # evaluate the polynomial coefficients
    x0 = xa[1]
    x1 = xa[2] - x0
    x2 = xa[3] - x0
    y0 = ya[1]
    y1 = ya[2] - y0
    y2 = ya[3] - y0
    t0 = 1.0/(x1*x2*(x2-x1))
    t1 = x1*y2
    t2 = x2*y1
    c1 = x2*t2 - x1*t1
    c2 = t1 - t2
    t1 = x - x0
    
    # evaluate the polynomial
    return y0 + t0*t1*(c1 + c2*t1)
end