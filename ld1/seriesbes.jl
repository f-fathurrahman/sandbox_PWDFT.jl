function seriesbes!(fun, r, r2, npt, xc)
#=
assume that the input function has the form
    xc(1) + xc(2)*r(n) + xc(3)*r(n)**2
and finds the two coefficients.
works with KKR3 beta functions
=#

#=
integer :: &
       npt,  &        ! the number of points  
       npt2          ! intermediate point

real(DP) :: &
       fun(npt),   &  ! the function
       r(npt),     &  ! the mesh     
       r2(npt),    &  ! the mesh     
       xc(4)         ! the coefficients
=#
    if npt < 3
        error("Need at least 3 points")
    end
    npt2 = floor(Int64, npt/2 + 1)
    # xc(1) = 0.5_dp*( fun(1) - xc(3)*r2(1) + fun(npt) - xc(3)*r2(npt) )

    xc[3] = ( ( fun[1] - fun[npt2] ) / ( r[1] - r[npt2] ) - 
              ( fun[npt] - fun[npt2] ) / ( r[npt] - r[npt2]) ) / ( r[1] - r[npt] )
    xc[1] = fun[1]
    xc[2] = ( fun[npt] - fun[npt2] ) / ( r[npt] - r[npt2] ) - xc[3] * ( r[npt] +  r[npt2] )
    xc[4] = 0.0
    return
end
