# !INPUT/OUTPUT PARAMETERS:
#   lmax : maximum angular momentum (in,integer)
#   zflm : coefficients of complex spherical harmonic expansion
#          (in,complex((lmax+1)**2)))
#   rflm : coefficients of real spherical harmonic expansion
#          (out,real((lmax+1)**2)))
function z_to_rf_lm!(lmax::Int64, zflm::Array{ComplexF64}, rflm::Array{Float64})

    # INTEGER, INTENT(in) :: lmax
    # COMPLEX(8), INTENT(in) :: zflm(*)
    # REAL(8), INTENT(out) :: rflm(*)

    c1 = 0.7071067811865475244  # 1/sqrt(2)
    lm1 = 0
    for l in 0:lmax
        lm2 = lm1 + 2*(l+1)
        for m in -l:-1
            lm1 = lm1+1
            lm2 = lm2-1
            if mod(m,2) != 0
                rflm[lm1] = -c1*( aimag(zflm[lm1]) + aimag(zflm[lm2]) )
            else
                rflm[lm1] =  c1*( aimag(zflm[lm2]) - aimag(zflm[lm1]) )
            end
        end
        lm1 = lm1 + 1
        lm2 = lm2 - 1
        rflm[lm1] = zflm[lm1]
        for m in 1:l
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            if mod(m,2) != 0
                rflm[lm1] = c1*(dble(zflm[lm1]) - dble(zflm[lm2]))
            else
                rflm[lm1] = c1*(dble(zflm[lm1]) + dble(zflm[lm2]))
            end
        end
    end
    return
end
