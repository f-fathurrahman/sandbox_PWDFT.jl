function series_besselj!(fun, r, r2, npt, xc)
#=
     assume that the input function has the form
        xc(1) + xc(2)*r(n) + xc(3)*r(n)**2
     and finds the two coefficients. works with KKR3 beta functions
=#

    @assert npt >= 3
    
    # intermediate point
    npt2 = floor(Int64, npt/2) + 1

    xc[3] = ( (fun[1] - fun[npt2]) / (r[1] - r[npt]) - 
              (fun[npt] - fun[npt2]) / (r[npt] - r[npt2]) ) / ( r[1] - r[npt] )
    xc[1] = fun[1]
    xc[2] = ( fun[npt] - fun[npt2] ) / ( r[npt] - r[npt2] ) - xc[3] * ( r[npt] + r[npt2] )
    xc[4] = 0.0
    
    return
end


