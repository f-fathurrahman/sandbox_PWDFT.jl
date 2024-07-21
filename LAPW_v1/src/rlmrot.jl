# Use boolean p?
function rlmrot!( p::Int64, ang::Vector{Float64}, lmax::Int64, D )
    #=
! !INPUT/OUTPUT PARAMETERS:
!   p    : if p=-1 then the rotation matrix is improper (in,integer)
!   ang  : Euler angles; alpha, beta, gamma (in,real(3))
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   d    : real spherical harmonic rotation matrix (out,real(ld,*))
! !DESCRIPTION:
!   Returns the rotation matrix in the basis of real spherical harmonics given
!   the three Euler angles, $(\alpha,\beta,\gamma)$, and the parity, $p$, of the
!   rotation.
=#

    SQTWO = sqrt(2)

    lmi = OffsetArray( zeros(Int64, 2*lmax+1), -lmax:lmax )
    ca = zeros(Float64, lmax)
    sa = zeros(Float64, lmax)
    cg = zeros(Float64, lmax)
    sg = zeros(Float64, lmax)
    
    ld = size(D, 1)
    dY = zeros(Float64, ld, ld)
    
    @assert lmax >= 0

    # generate the complex spherical harmonic rotation matrix about the y-axis
    ylmroty!( ang[2], lmax, dY)
    for m1 in 1:lmax
        ca[m1] = cos(m1*ang[1])
        sa[m1] = sin(m1*ang[1])
        cg[m1] = cos(m1*ang[3])
        sg[m1] = sin(m1*ang[3])
    end
    lm1 = 0
    for l in 0:lmax
        for m1 in -l:l
            lm1 = lm1 + 1
            lmi[m1] = lm1
        end
        lm0 = lmi[0]
        #println("lm0 = ", lm0)
        #@info "typeof(lm0) = $(typeof(lm0))"
        D[lm0,lm0] = dY[lm0,lm0]
        #@info "Pass here 46"
        for m1 in 1:l
            if mod(m1,2) == 0
                s1 = 1.0
            else
                s1 = -1.0
            end
            t1 = SQTWO * dY[lm0,lmi[m1]]
            t2 = SQTWO * dY[lmi[m1],lm0]
            D[lmi[m1], lm0] = s1*t1*ca[m1]
            D[lm0, lmi[m1]] = s1*t2*cg[m1]
            D[lmi[-m1], lm0] = -t1*sa[m1]
            D[lm0, lmi[-m1]] = t2*sg[m1]
            for m2 in 1:l
                if mod(m2,2) == 0
                    s2 = 1.0
                else
                    s2 = -1.0
                end
                t1 = ca[m1]*cg[m2]
                t2 = sa[m1]*sg[m2]
                t3 = sa[m1]*cg[m2]
                t4 = ca[m1]*sg[m2]
                t5 = dY[lmi[-m1], lmi[-m2]]
                t6 = s1*dY[lmi[m1], lmi[-m2]]
                t7 = t5 + t6
                t8 = t5 - t6
                D[lmi[ m1], lmi[ m2]] = s1*s2*(t1*t7 - t2*t8)
                D[lmi[ m1], lmi[-m2]] = s1*(t3*t8 + t4*t7)
                D[lmi[-m1], lmi[ m2]] = -s2*(t3*t7 + t4*t8)
                D[lmi[-m1], lmi[-m2]] = t1*t8 - t2*t7
            end
        end
    end

    # apply inversion if required
    if p == -1
        for l in range(1, lmax, step=2)
            lm1 = l^2 + 1
            lm2 = lm1 + 2*l
            D[lm1:lm2,lm1:lm2] = -D[lm1:lm2,lm1:lm2]
        end
    end

    return
end  # function
