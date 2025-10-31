# output are magmt and magir
function maginit!(
    atoms, rhomt, rhoir, ncmag, ndmag, bfcmt, bfieldc,
    magmt, magir
)

    # magnetization as fraction of density
    fmr = 0.15 # hard-coded
    v = zeros(Float64, 3)

    Natoms = atoms.Natoms
    # muffin-tin part
    for ia in 1:Natoms
        np = size(rhomt[ia], 1)
        v[:] = bfcmt[:,ia] + bfieldc[:]
        t1 = sqrt( v[1]^2 + v[2]^2 + v[3]^2 )
        if t1 > 1e-8
            t1 = -fmr/t1
            v[:] = t1*v[:]
            if !ncmag
                v[1] = v[3] # use z-component for collinear magn
            end
            for idm in 1:ndmag
                t1 = v[idm]
                magmt[ia][1:np] = t1*rhomt[ia][1:np]
            end
        else
            magmt[ia][1:np] = 0.0
        end
    end
    # interstitial magnetization
    v[:] = bfieldc[:]
    t1 = sqrt( v[1]^2 + v[2]^2 + v[3]^2 )
    if t1 > 1e-8
        t1 = -fmr/t1
        v[:] = t1*v[:]
        if !ncmag
            v[1] = v[3]
        end
        for idm in 1:ndmag
            t1 = v[idm]
            magir[:,idm] = t1*rhoir[:]
        end
    else
        magir[:,:] = 0.0
    end

    return
end