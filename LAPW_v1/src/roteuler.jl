
function roteuler!(R, angles)
#=INPUT/OUTPUT PARAMETERS:
!   rot : rotation matrix (in,real(3,3))
!   ang : Euler angles (alpha, beta, gamma) (out,real(3))
=#

    SMALL = 1e-8
    detR = det(R)
    if abs(detR - 1.0) > SMALL
        error("Matrix is improper or not unitary")
    end
    
    # No need for atan2 here.
    # In Julia atan function with two arguments (y,x) is atan2 in Fortran

    if (abs(R[3,1]) > SMALL) || (abs(R[3,2]) > SMALL)
        angles[1] = atan(R[3,2], R[3,1])
        if abs(R[3,1]) > abs(R[3,2])
            angles[2] = atan( R[3,1]/cos(angles[1]), R[3,3])
        else
            angles[2] = atan( R[3,2]/sin(angles[1]), R[3,3])
        end
        angles[3] = atan( R[2,3], -R[1,3])
    else
        angles[1] = atan(R[1,2], R[1,1])
        if R[3,3] > 0.0
            angles[2] = 0.0
            angles[3] = 0.0
        else
            angles[2] = pi
            angles[3] = pi
        end
    end
    return
end

# ccall( (:roteuler_, LIBLAPW), Nothing, (Ptr{Float64}, Ptr{Float64}), R, ang)
# Multiply with -1 if the determinant is -1 to make it proper rotation.