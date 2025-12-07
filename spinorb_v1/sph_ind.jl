#=
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
=#

function sph_ind( l::Int64, j::Int64, m::Int64, spin::Int64 )
    #=
    This function calculates the m index of the spherical harmonic
    in a spinor with orbital angular momentum l, total angular 
    momentum j, projection along z of the total angular momentum m+-1/2. 
    Spin selects the up (spin=1) or down (spin=2) coefficient.
    =#
    @assert spin in [1,2]
    @assert m >= (-l - 1)
    @assert m <= l
    #
    if abs(j - l - 0.5) < 1e-8
        if spin == 1
            return m
        end
        if spin == 2
            return m + 1
        end
    elseif abs(j - l + 0.5) < 1e-8 # check if j - l == 1/2
        if m < (-l + 1)
            return 0
        else
            if spin == 1
                return m - 1
            end
            if spin == 2
                return m
            end
        end
    else
        error("l=$l and j=$j are not compatible")
    end
    #
    if (sph_ind < -l) || (sph_ind > l)
        return 0
    end
end

