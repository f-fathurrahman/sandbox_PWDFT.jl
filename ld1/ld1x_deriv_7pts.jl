function deriv_7pts(f, ik, rc, h)
    #=
! ---------------------------------------------------------------
  !      evaluates the first derivative of function f, the function
  !      is given numerically on logarithmic mesh r.
  !      nm = dimension of mesh
  !      ik = integer : position of the point in which the derivative
  !      will be evaluated.
  !      h is distance between x(i) and x(i+1) where
  !      r(i) = exp(x(i))/znesh & r(j) = exp(x(j))/znesh
  !
  use kinds, only : DP
  implicit none
  integer :: a(7),n,ik,i
  real(DP) :: f(ik+3),rc,h,sum,deriv_7pts
  !      coefficients for the formula in abramowitz & stegun p.914
  data a/-12,108,-540,0,540,-108,12/
    =#

    a = [-12, 108, -540, 0, 540, -108, 12]

    # formula for linear mesh 
    ss = 0
    for i in 1:7
        ss += a[i]*f[i-4+ik]
    end
    res = sum/(720.0*h)
    # transform to logarithmic mesh
    res = res/rc
    return res
end

