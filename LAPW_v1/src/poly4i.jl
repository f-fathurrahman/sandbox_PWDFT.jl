function poly4i(xa, ya, x)
    # arguments
    #real(8), intent(in) :: xa(4),ya(4),x

    # evaluate the polynomial coefficients
    x0 = xa[1]
    x1 = xa[2] - x0
    x2 = xa[3] - x0
    x3 = xa[4] - x0
    y0 = ya[1]
    y1 = ya[2] - y0
    y2 = ya[3] - y0
    y3 = ya[4] - y0
    t4 = x1 - x2
    t5 = x1 - x3
    t6 = x2 - x3
    t1 = x1*x2*y3
    t2 = x2*x3*y1
    t3 = x1*x3
    t0 = 1.0/(x2*t3*t4*t5*t6)
    t3 = t3*y2
    c3 = t1*t4 + t2*t6 - t3*t5
    t4 = x1^2
    t5 = x2^2
    t6 = x3^2
    c2 = t1*(t5 - t4) + t2*(t6 - t5) + t3*(t4 - t6)
    c1 = t1*(x2*t4 - x1*t5) + t2*(x3*t5 - x2*t6) + t3*(x1*t6 - x3*t4)
    t1 = x - x0
    # integrate the polynomial
    return t1*(y0 + t0*t1*(0.5*c1 + t1*(c2/3 + 0.25*c3*t1) ) )
end