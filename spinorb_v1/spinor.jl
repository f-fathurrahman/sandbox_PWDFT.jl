#=
!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
=#
function spinor( l, j, m, spin )
    #=
    This function calculates the numerical coefficient of a spinor
    with orbital angular momentum l, total angular momentum j, 
    projection along z of the total angular momentum m+-1/2. Spin selects
    the up (spin=1) or down (spin=2) coefficient.
    =#

    #=
    l:    orbital angular momentum
    m:    projection of the total angular momentum+-1/2
    spin: 1 or 2 select the component
    j:    total angular momentum
    =#

    @assert spin in [1,2]
    @assert m >= (-l - 1)
    @assert m <= l
    
    denom = 1.0 / (2*l + 1)
    
    if abs(j - l - 0.5) < 1e-8
        #
        if spin == 1
            return sqrt((l + m + 1)*denom)
        end
        #
        if spin == 2
            return sqrt((l-m)*denom)
        end
        #
    elseif abs(j - l + 0.5) < 1e-8
        #
        if m < (-l + 1)
            return 0.0
        else
            if spin == 1
                return sqrt((l - m + 1)*denom)
            end
            if spin == 2
                return -sqrt((l + m)*denom)
            end
        end
        #
    else
        error("j and l not compatible")
    end

    error("Should not go here")
end
