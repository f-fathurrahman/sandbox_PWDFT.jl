#=
! !INPUT/OUTPUT PARAMETERS:
!   m    : order of derivatve (in,integer)
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   djl  : array of RETURN ed values (out,real(0:lmax))
=#
function sbesseldm!(m, lmax, x, djl)
    #=
    IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: m,lmax
  REAL(8), intent(in) :: x
  REAL(8), intent(out) :: djl(0:lmax)
  ! local variables
  INTEGER i,j,l,i0
  REAL(8) t1,sum,x2
  INTEGER a(0:maxm+1),a1(0:maxm+1)
  INTEGER b(0:maxm+1),b1(0:maxm+1)
  ! automatic arrays
  REAL(8) jl(0:lmax+1)
  ! external functions
  REAL(8) factnm,factr
  external factnm,factr
    =#

    MAXM = 6
    MAXNS = 20

    if (m < 0) || (m > MAXM)
        error("m out of range")
    end

    if (lmax < 0) || (lmax > 39)
        error("lmax out of range")
    end
    
    if (x < 0.0) || (x > 1e5)
        error("x out of range")
    end
  
    if m == 0
        sbessel!(lmax, x, djl)
        return
    end

    a = OffsetArray( zeros(Int64, MAXM+2), 0:MAXM+1 )
    a1 = OffsetArray( zeros(Int64, MAXM+2), 0:MAXM+1 )
    b = OffsetArray( zeros(Int64, MAXM+2), 0:MAXM+1 )
    b1 = OffsetArray( zeros(Int64, MAXM+2), 0:MAXM+1 )
    jl = OffsetArray( zeros(Float64, lmax+2), 0:lmax+1 )

    if x > 1.0
        sbessel!(lmax+1, x, jl)
        for l in 0:lmax
            a[1:m+1] .= 0
            a[0] = 1
            a1[0:m+1] .= 0
            for i in 1:m
                b[0] = 0
                b1[0] = 0
                for j in 0:i
                    b[j+1] = a[j]*(l - j)
                    b1[j+1] = -a1[j]*(j + l + 2)
                end
                for j in 0:i
                    b1[j] = b1[j] - a[j]
                    b[j] = b[j] + a1[j]
                end 
                a[0:i+1] .= b[0:i+1]
                a1[0:i+1] .= b1[0:i+1]
            end
            t1 = 1.0
            ss = a[0] * jl[l] + a1[0] * jl[l+1]
            for i in 1:(m+1)
                t1 *= x # update
                ss += ( a[i] * jl[l] + a1[i] * jl[l+1] ) / t1 # update
            end
            djl[l] = ss
        end 
    else
        x2 = x^2
        for l in 0:lmax
            i0 = max( floor(Int64, (m-l+1)/2), 0) # use round?
            j = 2*i0 + l - m
            if j == 0
                t1 = 1.0
            else
                t1 = x^j
            end 
            t1 = factr(j+m,j)*t1/(factnm(i0,1)*factnm(j+l+m+1,2) *  (-2)^i0)
            ss = t1
            for i in (i0+1):MAXNS
                j = 2*i + l
                t1 = -t1 * (j-1)*j * x2 / ( (j-l)*(j-m-1)*(j-m)*(j+l+1) )
                if abs(t1) <= 1.e-40
                    @goto LABEL10 # break
                end
                ss += t1 # update
            end
            @label LABEL10 # continue
            djl[l] = ss
        end 
    end # if
    
    return
end
