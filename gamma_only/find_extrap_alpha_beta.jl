function find_extrap_alpha_beta( tau, tauold )
    #
    # ... This routine finds the best coefficients alpha0 and beta0 so that
    #
    # ...    | tau(t+dt) - tau' | is minimum, where
    #
    # ...    tau' = tau(t) + alpha0 * ( tau(t) - tau(t-dt) )
    # ...                  + beta0 * ( tau(t-dt) -tau(t-2*dt) )
    #
    # REAL(DP) :: alpha0, beta0, tau(3,nat), tauold(3,nat,3)
    # REAL(DP) :: a11, a12, a21, a22, b1, b2, c, det
    
    # IF ( history <= 2 ) RETURN
    #
    # ... solution of the linear system
    #
    a11 = 0.0
    a12 = 0.0
    a21 = 0.0
    a22 = 0.0
    b1  = 0.0
    b2  = 0.0
    c   = 0.0
    
    Natoms = size(tau,2)
    for ia in 1:Natoms
        for ipol in 1:3
            a11 = a11 + ( tauold[ipol,ia,1] - tauold[ipol,ia,2] )^2
            a12 = a12 + ( tauold[ipol,ia,1] - tauold[ipol,ia,2] ) * 
                        ( tauold[ipol,ia,2] - tauold[ipol,ia,3] )
            a22 = a22 + ( tauold[ipol,ia,2] - tauold[ipol,ia,3] )^2
            b1 = b1 - ( tauold[ipol,ia,1] - tau[ipol,ia] ) * 
                      ( tauold[ipol,ia,1] - tauold[ipol,ia,2] )
            b2 = b2 - ( tauold[ipol,ia,1] - tau[ipol,ia] ) * 
                      ( tauold[ipol,ia,2] - tauold[ipol,ia,3] )
            c = c + ( tauold[ipol,ia,1] - tau[ipol,ia] )^2
        end
    end

    SMALL = 1e-15

    a21 = a12
    D = a11 * a22 - a12 * a21
    if D < -SMALL
       alpha0 = 0.0
       beta0  = 0.0
       println("SMALL determinant: ", D)
    end
    
    #
    # ... case D > 0:  a well defined minimum exists
    #
    if D > SMALL
        alpha0 = ( b1 * a22 - b2 * a12 ) / D
        beta0  = ( a11 * b2 - a21 * b1 ) / D
    else
        #
        # ... case det = 0 : the two increments are linearly dependent,
        # ...                chose solution with alpha = b1 / a11 and beta = 0
        # ...                ( discard oldest configuration )
        #
        println("SMALL D v2 = ", D)
        alpha0 = 0.0
        beta0  = 0.0
        if a11 != 0.0
            alpha0 = b1 / a11
        end
    end

  return alpha0, beta0

end