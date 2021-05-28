# !INPUT/OUTPUT PARAMETERS:
#   nr     : number of radial mesh points (in,integer)
#   nri    : number of points on inner part of muffin-tin (in,integer)
#   ld     : leading dimension (in,integer) REMOVED
#   rl     : r^l on the radial mesh (in,real(ld,-lmaxo-1:lmaxo+2))
#   wpr    : weights for partial integration on radial mesh (in,real(4,nr))
#   zrhomt : muffin-tin charge density (in,complex(*))
#   zvclmt : muffin-tin Coulomb potential (out,complex(*))
function zpotclmt!( mt_vars, nr, nri, rl, wpr, zrhomt, zvclmt )
    # USE m_constants, ONLY: fourpi
    # USE m_muffin_tins, ONLY: lmaxo, lmmaxo, lmmaxi, lmaxi
    # IMPLICIT NONE 
    # ! arguments
    # INTEGER, intent(in) :: nr,nri,ld
    # REAL(8), intent(in) :: rl(ld,-lmaxo-1:lmaxo+2),wpr(4,nr)
    # COMPLEX(8), intent(in) :: zrhomt(*)
    # COMPLEX(8), intent(out) :: zvclmt(*)
    # ! local variables
    # INTEGER :: nro,iro,ir
    # INTEGER :: l,l1,l2,l3
    # INTEGER :: m,lm,npi,i
    # REAL(8) :: r1,r2,t0,t1,t2,t3,t4
    # ! automatic arrays
    # REAL(8) :: f1(nr),f2(nr),f3(nr),f4(nr),f5(nr)

    lmaxo = mt_vars.lmaxo
    lmaxi = mt_vars.lmaxi
    lmmaxo = mt_vars.lmmaxo
    lmmaxi = mt_vars.lmmaxi

    nro = nr-nri
    iro = nri+1
    npi = lmmaxi*nri
    lm = 0
    for l in 0:lmaxi
        l1 = l + 2
        l2 = -l + 1
        l3 = -l - 1
        t0 = 4*pi/(2*l + 1)
        for m in -l:l
            lm = lm+1
            i = lm
            for ir in 1:nri
                t1 = real(zrhomt[i])
                t2 = imag(zrhomt[i])
                r1 = rl[ir,l1]
                r2 = rl[ir,l2]
                f1[ir] = t1*r1
                f2[ir] = t2*r1
                f3[ir] = t1*r2
                f4[ir] = t2*r2
                i = i + lmmaxi
            end
            for ir in iro:nr
                t1 = real(zrhomt[i])
                t2 = imag(zrhomt[i])
                r1 = rl[ir,l1]
                r2 = rl[ir,l2]
                f1[ir] = t1*r1
                f2[ir] = t2*r1
                f3[ir] = t1*r2
                f4[ir] = t2*r2
                i = i + lmmaxo
            end
            splintwp!(nr, wpr, f1, f5)
            splintwp!(nr, wpr, f2, f1)
            splintwp!(nr, wpr, f3, f2)
            splintwp!(nr, wpr, f4, f3)
            t1 = f2[nr]
            t2 = f3[nr]
            i = lm
            for ir in 1:nri
                r1 = t0*rl[ir,l3]
                r2 = t0*rl[ir,l]
                t3 = r1*f5[ir] + r2*(t1 - f2[ir])
                t4 = r1*f1[ir] + r2*(t2 - f3[ir])
                zvclmt[i] = t3 + im*t4
                i = i + lmmaxi
            end
            for ir in iro:nr
                r1 = t0*rl[ir,l3]
                r2 = t0*rl[ir,l]
                t3 = r1*f5[ir] + r2*(t1 - f2[ir])
                t4 = r1*f1[ir] + r2*(t2 - f3[ir])
                zvclmt[i] = t3 + im*t4
                i = i + lmmaxo
            end 
        end
    end
    
    for l in lmaxi+1:lmaxo
        l1 = l + 2
        l2 = -l + 1
        l3 = -l - 1
        t0 = 4*pi/(2*l+1)
        for m in -l:l
            lm = lm + 1
            i = npi + lm
            for ir in iro:nr
                t1 = real(zrhomt[i])
                t2 = imag(zrhomt[i])
                r1 = rl[ir,l1]
                r2 = rl[ir,l2]
                f1[ir] = t1*r1
                f2[ir] = t2*r1
                f3[ir] = t1*r2
                f4[ir] = t2*r2
                i = i + lmmaxo
            end
            idx = iro:iro+nro-1
            @views splintwp( nro, wpr[:,idx], f1[idx], f5[idx] )
            @views splintwp( nro, wpr[:,idx], f2[idx], f1[idx] )
            @views splintwp( nro, wpr[:,idx], f3[idx], f2[idx] )
            @views splintwp( nro, wpr[:,idx], f4[idx], f3[idx] )
            t1 = f2[nr]
            t2 = f3[nr]
            i = npi + lm
            for ir in iro:nr
                r1 = t0*rl[ir,l3]
                r2 = t0*rl[ir,l]
                t3 = r1*f5[ir] + r2*(t1 - f2[ir])
                t4 = r1*f1[ir] + r2*(t2 - f3[ir])
                zvclmt[i] = t3 + im*t4
                i = i + lmmaxo
            end
        end
    end
    return
end