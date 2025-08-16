#
# Copyright (C) 2001-2007 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

# no of points are taken to be from r and not from jl
function qe_sph_bes!( l::Int64, q::Float64, r, jl )

    eps14 = 1e-14

    # case q=0
    if abs(q) < eps14
        if l == -1
            error("Invalid input: l=$l")
        elseif l == 0
            fill!(jl, 1.0)
        else
            fill!(jl, 0.0)
        end # if
        return
    end #if 

    # case l = -1

    if l == - 1
        if abs(q * r[1]) < eps14
            error("Invalid inputL l=$(l) q=$(q)")
        end
        @views jl[:] = cos(q*r[:]) / (q*r[:])
        return
    end # if

    # series expansion for small values of the argument
    # ir0 is the first grid point for which q*r(ir0) > xseries
    # notice that for small q it may happen that q*r(Npoints) < xseries !
    Npoints = length(r)
    XSERIES = 0.05

    ir0 = Npoints + 1
    for ir in 1:Npoints
        if abs(q*r[ir]) > XSERIES
            ir0 = ir
            break
        end
    end

    for ir in 1:(ir0-1)
        x = q * r[ir]
        if l == 0
            xl = 1.0
        else
            xl = x^l
        end # if
        jl[ir] = xl/besselj_semifact(2*l+1) * ( 1.0 - x^2/1.0/2.0/(2.0*l+3) *
                                      ( 1.0 - x^2/2.0/2.0/(2.0*l+5) *
                                      ( 1.0 - x^2/3.0/2.0/(2.0*l+7) *
                                      ( 1.0 - x^2/4.0/2.0/(2.0*l+9) ) ) ) )
    end # do

    if ir0 > Npoints
        return
    end

    idxs = ir0:Npoints
    if l == 0
        #
        @. jl[idxs] = sin(q*r[idxs]) / (q*r[idxs])
        #
    elseif l == 1
        @. jl[idxs] = ( sin(q*r[idxs]) / (q*r[idxs] ) - cos(q*r[idxs] ) ) / (q*r[idxs] )
        #
    elseif l == 2

        @. jl[idxs] = ( (3/(q*r[idxs]) - (q*r[idxs])) * sin(q*r[idxs]) - 3 * cos(q*r[idxs]) ) / (q*r[idxs])^2

    elseif l == 3

        @. jl[idxs] = (sin(q*r[idxs]) * (15 / (q*r[idxs]) - 6*(q*r[idxs]) ) + 
                      cos(q*r[idxs]) * ( (q*r[idxs])^2 - 15) ) / (q*r[idxs])^3

    elseif l == 4

        @. jl[idxs] = (sin(q*r[idxs]) * (105 - 45 * (q*r[idxs])^2 +
                     (q*r[idxs])^4) + cos(q*r[idxs]) * (10*(q*r[idxs])^3 - 105*(q*r[idxs])) ) / (q*r[idxs])^5
  
    elseif l == 5

        @. jl[idxs] = (-cos(q*r[idxs]) - 
                  (945*cos(q*r[idxs])) / (q*r[idxs])^4 +
                  (105*cos(q*r[idxs])) / (q*r[idxs])^2 +
                  (945*sin(q*r[idxs])) / (q*r[idxs])^5 -
                  (420*sin(q*r[idxs])) / (q*r[idxs])^3 +
                  ( 15*sin(q*r[idxs])) / (q*r[idxs]) ) / (q*r[idxs])

    elseif l == 6

        @. jl[idxs] = ((-10395*cos(q*r[idxs])) / (q*r[idxs])^5 +
                  (  1260*cos(q*r[idxs])) / (q*r[idxs])^3 -
                  (    21*cos(q*r[idxs])) / (q*r[idxs]) - sin(q*r[idxs]) +
                  ( 10395*sin(q*r[idxs])) / (q*r[idxs])^6 -
                  (  4725*sin(q*r[idxs])) / (q*r[idxs])^4 +
                  (   210*sin(q*r[idxs])) / (q*r[idxs])^2 ) / (q*r[idxs])
    else

        error("Invalid input: l=$(l)")

    end

    return

end # function


function besselj_semifact(n)
    # semifact(n) = n!!   
    res = 1
    for i in range(n, stop=1, step=-2)
        res = i*res
    end
    return res
end
