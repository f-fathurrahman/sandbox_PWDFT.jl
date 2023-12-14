#
# Copyright (C) 2009 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

# assumption: points in r are ordered from small to large
function _rgrad_find_right_point(i, r, Δ)
    Nr = size(r, 1)
    jidx = -1
    is_found = false
    for j in (i+1):Nr
        if r[j] > ( r[i] + Δ )
            jidx = j
            is_found = true
            break
        end
    end
    return is_found, jidx
end


function _rgrad_find_left_point(i, r, Δ)
    Nr = size(r, 1)
    kidx = -1
    is_found = false
    for k in range(i-1, stop=1, step=-1)
        if r[k] < (r[i] - Δ)
            kidx = k
            is_found = true
        end
    end
    return is_found, kidx
end

function radial_gradient_coarse!(r, f, gf)
    
    Nr = size(r, 1)
    Δ = 1e-5
    imin = 1
    for i in 2:Nr # LABEL points
        is_found_right, j = _rgrad_find_right_point(i, r, Δ)
        is_found_left, k = _rgrad_find_left_point(i, r, Δ)
        if is_found_right && is_found_left
            num1 = (r[j] - r[i])^2 * (f[k] - f[i]) - (r[k] - r[i])^2*( f[j] - f[i] )
            denum1 = (r[j] - r[i]) * (r[k] - r[i]) * (r[j] - r[k])
            gf[i] = num1/denum1
        end
        if !is_found_left
            imin = i
        end
        # set gf[i] to zero for the moment
        # It will be calculated by interpolation below
        if !is_found_right
            gf[i] = 0.0
        end
    end # do points
    #println("Sum gf after loop = ", sum(gf))
    #println("imin = ", imin)

    #
    # In the first imin points the previous formula cannot be
    # used. We interpolate with a polynomial the points already found
    # and extrapolate in the points from 1 to imin.
    # Presently we fit 5 points with a 3rd degree polynomial.
    #
    npoint = 5
    b = zeros(Float64, npoint)
    faux = zeros(Float64, npoint+1)
    raux = zeros(Float64, npoint+1)
    faux[1] = gf[imin+1]
    raux[1] = r[imin+1]
    j = imin + 1
    for k in 2:npoint # LABEL points fit
        for i in j:(Nr-1)
            if r[i] > (r[imin+1] + (k-1)*Δ)
                faux[k] = gf[i]
                raux[k] = r[i]
                j = i + 1
                # break current loop, next iteration for points fit
                break
            end #if
        end #do
    end #do points_fit

    #println("Before fit_pol! b = ", b)    
    _rgrad_fit_pol!(raux, faux, npoint, 3, b)
    #println("After fit_pol! b = ", b)

    # evaluate the polynomial
    for i in 1:imin
        gf[i] = b[1] + r[i]*(b[2] + r[i]*(b[3] + r[i]*b[4]) )
    end
    return
end

# This routine finds the coefficients of the least-square polynomial which 
# interpolates the n input data points.
function _rgrad_fit_pol!( xdata, ydata, n, degree, b)
    bmat = zeros(Float64, degree+1, degree+1)
    amat = zeros(Float64, degree+1, n)
    amat[1,:] .= 1.0
    for i in 2:(degree+1)
        for j in 1:n
            amat[i,j] = amat[i-1,j]*xdata[j]
        end
    end
    #
    for i in 1:(degree+1)
        b[i] = 0.0
        for k in 1:n
            b[i] += ydata[k]*xdata[k]^(i-1)
        end
    end
    #
    for i in 1:(degree+1)
        for j in 1:(degree+1)
            bmat[i,j] = 0.0
            for k in 1:n
                bmat[i,j] += amat[i,k]*amat[j,k]
            end
        end
    end
    
    #println("size bmat = ", size(bmat))
    #println("size b = ", size(b))
    # Solve the linear equation  
    @views LAPACK.gesv!(bmat, b[1:(degree+1)])  
    return
end


