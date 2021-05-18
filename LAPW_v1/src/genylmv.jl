# !INPUT/OUTPUT PARAMETERS:
#   lmax : maximum angular momentum (in,integer)
#   v    : input vector (in,real(3))
#   ylm  : array of spherical harmonics (out,complex((lmax+1)^2))
function genylmv!(lmax::Int64, v, ylm)
    # IMPLICIT NONE 
    # ! arguments
    # INTEGER, intent(in) :: lmax
    # REAL(8), intent(in) :: v(3)
    # COMPLEX(8), intent(out) :: ylm(*)
    # ! local variables
    # REAL(8), parameter :: eps=1.d-14
    # COMPLEX(8) z1
  
    @assert lmax >= 0
    @assert lmax <= 50

    SMALL = 1e-14

    ylm[1] = 0.28209479177387814347
    if lmax == 0
        return
    end

    r = sqrt( v[1]^2 + v[2]^2 + v[3]^2)
    st = 0.0; ct = 0.0; sp = 0.0; cp1 = 0.0 # initialize
    if r > SMALL
        t1 = v[3]/r
        if t1 >= 1.0
            st = 0.0
            ct = 1.0
        elseif t1 <= -1.0
            st = 0.0
            ct = -1.0
        else
            st = sqrt(1.0 - t1^2)
            ct = t1
        end
        if (abs(v[1]) > SMALL) || (abs(v[2]) > SMALL)
          t1 = 1.0/sqrt(v[1]^2 + v[2]^2)
          sp = t1*v[2]
          cp1 = t1*v[1]
        else
          sp = 0.0
          cp1 = 1.0
        end
    else 
        st = 0.0
        ct = 1.0
        sp = 0.0
        cp1 = 1.0
    end

    z1 = cp1 + im*sp
    ylm[3] = 0.48860251190291992159*ct
    ylm[4] = -0.34549414947133547927*st*z1
    ylm[2] = -conj(ylm[4])
  
    for l in 2:lmax
        lm1 = (l+1)^2
        lm2 = l^2
        lm3 = (l-1)^2+1
        lm4 = lm2+1
        ylm[lm1] = -st*sqrt((2*l+1)/(2*l)) * z1 * ylm[lm2]
        #
        if mod(l,2) == 0
          ylm[lm4] = conj(ylm[lm1])
        else
          ylm[lm4] = -conj(ylm[lm1])
        end
        #
        lm1 = lm1 - 1
        ylm[lm1] = ct*sqrt(2*l+1)*ylm[lm2]
        lm4 = lm4 + 1
        if mod(l-1,2) == 0
            ylm[lm4] = conj(ylm[lm1])
        else
            ylm[lm4] = -conj(ylm[lm1])
        end
        #
        t1 = ct*sqrt((2*l-1)*(2*l+1))
        t2 = sqrt((2*l+1)/(2*l-3))
        for m in range(l-2, stop=1, step=-1)
            lm1 = lm1 - 1
            lm2 = lm2 - 1
            lm3 = lm3 - 1
            lm4 = lm4 + 1
            t3 = 1.0/sqrt((l-m)*(l+m))
            t4 = t2*sqrt((l-m-1)*(l+m-1))
            ylm[lm1] = t3*(t1*ylm[lm2] - t4*ylm[lm3])
            if mod(m,2) == 0
              ylm[lm4] = conj(ylm[lm1])
            else
              ylm[lm4] = -conj(ylm[lm1])
            end
        end
        lm1 = lm1-1; lm2 = lm2-1; lm3 = lm3-1
        t3 = 1.0/l
        t4 = t2*(l-1)
        ylm[lm1] = t3*( t1*ylm[lm2] - t4*ylm[lm3] )
    end
    return
end
