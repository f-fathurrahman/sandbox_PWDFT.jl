#
# Copyright (C) 2010 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

# determines the wave-function in the first two points by
# series developement. It receives as input:
# lam the angular momentum 
# e   the energy in Ry (original in QE)
# b(0:3) the coefficients of a polynomial that interpolates the 
#        potential in the first three points
# grid  the mesh
# ze2   the zed of the mesh, changed to Zval
# in output solution(1:2) contains the solution in the first two points
function start_scheq!(l::Int64, e, b, grid::RadialGrid, Zval, solution)
    #
    #  set up constants and initialize
    #
    l1 = l + 1
    xl1 = l + 1.0
    x4l6 = 4.0*l + 6.0
    x6l12 = 6.0*l + 12.0
    x8l20 = 8.0*l + 20.0

    b0e = b[1] - e
    c1 = Zval/xl1
    c2 = (c1*Zval + b0e)/x4l6
    c3 = (c2*Zval + c1*b0e + b[2])/x6l12
    c4 = (c3*Zval + c2*b0e + c1*b[2] + b[3])/x8l20
    r = grid.r

    rr1 = ( 1.0 + r[1]*(c1 + r[1]*(c2 + r[1]*(c3 + r[1]*c4))) ) * r[1]^l1
    rr2 = ( 1.0 + r[2]*(c1 + r[2]*(c2 + r[2]*(c3 + r[2]*c4))) ) * r[2]^l1
    
    solution[1] = rr1/grid.sqrtr[1]
    solution[2] = rr2/grid.sqrtr[2]

    return
end
