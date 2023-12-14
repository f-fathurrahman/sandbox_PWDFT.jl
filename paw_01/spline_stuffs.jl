#
# Copyright (C) 2004-2006 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

function init_spline!(
    xdata, ydata,
    startu::Float64, startd::Float64,
    d2y
)
    ydim = size(ydata, 1)
    u = zeros(Float64, ydim)
    u[1] = startu
    d2y[1] = startd
    
    for i in 2:(ydim-1)
        #
        sig = ( xdata[i] - xdata[i-1] ) / ( xdata[i+1] - xdata[i-1] ) 
        p = sig * d2y[i- 1] + 2.0 
        d2y[i] = ( sig - 1.0 ) / p 
        u[i] = ( 6.0 * ( ( ydata[i+1] - ydata[i] ) / ( xdata[i+1] - xdata[i] ) - 
                         ( ydata[i] - ydata[i-1] ) / ( xdata[i] - xdata[i-1] ) ) /
                         ( xdata[i+1] - xdata[i-1] ) - sig * u[i-1] ) / p 
    end
    # TODO: make sympy expression of the above operation
    #
    d2y[ydim] = 0.0  
    for k in range(ydim-1, stop=1, step=-1)
        d2y[k] = d2y[k] * d2y[k+1] + u[k] 
    end
    #
    return
end

function spline_interp( xdata, ydata, d2y, x )     
    xdim = size(xdata, 1)
    klo = 1
    khi = xdim
    iloc = _spline_point_locate(xdata, x)
    println("iloc = ", iloc)
    klo = max(min(iloc,xdim - 1), 1)
    println("klo = ", klo)
    khi = klo + 1
    h = xdata[khi] - xdata[klo]
    a = ( xdata[khi] - x ) / h
    b = ( x - xdata[klo] ) / h
    res = a * ydata[klo] + b * ydata[khi] + ( ( a^3 - a ) * d2y[klo] + ( b^3 - b ) * d2y[khi] ) * h^2 / 6.0
    return res
end

function _spline_point_locate( xx::Vector{Float64}, x::Float64 )
    
    n = size(xx, 1)
    ascnd = xx[n] >= xx[1]
    #println("ascnd = ", ascnd)
    #
    jl = 0
    ju = n + 1

    # bisection?
    while true # main loop
        if (ju - jl) <= 1
            break
        end
        jm = floor(Int64, (ju + jl)/2)
        if ascnd == (x >= xx[jm])
            jl = jm
        else
            ju = jm
        end
    end
    # FIXME case jl or ju falls outside valid range?
    
    SMALL = eps(Float64)
    idx_return = -1 # set to an invalid value
    if abs(x - xx[1]) <= SMALL
        #
        idx_return = 1
        #
    elseif abs(x - xx[n]) <= SMALL
        #
        idx_return = n - 1
        #
    else 
        #
        idx_return = jl
    end

    return idx_return
end

