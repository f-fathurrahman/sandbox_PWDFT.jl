function calc_mag_moment!( atoms, pw, mt_vars, cfunir, elec_chgst )
#=
IMPLICIT NONE 
! local variables
INTEGER idm,is,ias,nr,nri
REAL(8) t1
! automatic arrays
REAL(8) fr(nrmtmax)
=#
    fr = zeros(Float64, nrmtmax)

    # find the muffin-tin moments
    fill!(mommttot, 0.0)
    for idm in 1:ndmag
        for ia in 1:Natoms
            isp = atm2species[ia]
            nr = nrmt[isp]
            nri = nrmti[isp]
            # extract the l=m=0 component from the muffin-tin magnetisation
            rf_mt_lm!(1, mt_vars, isp, magmt[ia][:,idm], fr)
            # integrate to the muffin-tin radius
            t1 = dot( wrmt[isp], fr[1:nr])
            mommt[ia][idm] = 4π * y00 * t1
            mommttot[idm] = mommttot[idm] + mommt[ia][idm]
        end
    end
    # find the interstitial moments
    for idm in 1:ndmag
        t1 = dot( magir[:,idm], cfunir )
        momir[idm] = t1*CellVolume/Npoints
    end
    @. momtot[:] = mommttot[:] + momir[:]
    # total moment magnitude
    if ncmag
        momtotm = sqrt(momtot[1]^2 + momtot[2]^2 + momtot[3]^2)
    else 
        momtotm = abs(momtot[1])
    end
    return
end