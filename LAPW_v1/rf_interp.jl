# !INPUT/OUTPUT PARAMETERS:
#   ni : number of input points (in,integer)
#   xi : input abscissa array (in,real(ni))
#   fi : input data array (in,real(ni)
#   no : number of output points (in,integer)
#   xo : output abscissa array (in,real(ni))
#   fo : output interpolated function (out,real(no))
function rf_interp!(ni, xi, fi, no, xo, fo)
  
    @assert ni > 0
    @assert no > 0
  
    # HUH? Probably not needed
    if ni==1
        fo[:] .= fi[1]
        return
    end

    # compute the spline coefficients
    cf = zeros(Float64,3,ni)
    spline!(ni, xi, fi, cf)
    
    # evaluate spline at output points
    for l in 1:no
        x = xo[l]
        i = 1
        # XXX: Need a better algoritm
        while true
            is_in_interval = (x >= xi[i]) & (x <= xi[i+1])
            if is_in_interval
                dx = x - xi[i]
                fo[l] = fi[i] + dx*( cf[1,i] + dx*( cf[2,i] + dx*cf[3,i] ) )
                break
            else
                i = i + 1
            end
            if i == ni
                error("Outside the interval")
                break
            end
        end
    end
    return
end