function deriv2_7pts(f, ik, rc, h)
#=
evaluates the second derivative of function f, the function
is given numerically on logarithmic mesh r.

nm = dimension of mesh
ik = integer : position of the point in which the derivative
will be evaluated.

h is distance between x(i) and x(i+1) where
r(i) = exp(x(i))/znesh & r(j) = exp(x(j))/znesh

  integer :: a(7),n,nm,i,ik
  real(DP) :: f(ik+3),
  rc,h,sum,sum1,deriv_7pts,deriv2_7pts
  !      coefficients for the formula in abramowitz & stegun p.914
  !      the formula is used for 7 points.
  data a/4,-54,540,-980,540,-54,4/   ! these are coefficients
=#

    a = [4, -54, 540, -980, 540, -54, 4]

    # formula for linear mesh 
    ss = 0.0
    for i in 1:7
        ss = ss + a[i] * f[i-4+ik]
    end
    ss = 2.0*ss/(720.0 * h^2)
    # transform to logarithmic mesh
    ss1 = deriv_7pts(f, ik, rc, h)
    return ss/(rc*rc) - ss1/rc
end