function mtdmin( atoms::Atoms, rmt::Array{Float64,1} )
    v1 = zeros(Float64,3)
    v2 = zeros(Float64,3)
    is = 1
    js = 1
    dmin = 1.0e6

    LatVecs = atoms.LatVecs
    epslat = 1e-6
    atpos = atoms.positions
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    
    for i1 in -1:1, i2 in -1:1, i3 in -1:1
        #
        @views v1[:] = i1*LatVecs[:,1] + i2*LatVecs[:,2] + i3*LatVecs[:,3]
        #
        for ka in 1:Natoms
            ks = atm2species[ka]
            #
            @views v2[:] = v1[:] + atpos[:,ka]
            #
            for la in 1:Natoms
                #
                ls = atm2species[la]
                t1 = rmt[ks] + rmt[ls]
                if ( (i1 != 0) || (i2 != 0) || (i3 != 0) || (ks != ls) ||  (ka != la) )
                    t2 = sqrt( (v2[1] - atpos[1,la])^2 +
                               (v2[2] - atpos[2,la])^2 +
                               (v2[3] - atpos[3,la])^2 )
                    t3 = t2 - t1
                    if t3 < (dmin - epslat)
                        is = ks
                        js = ls
                        dmin = t3
                    end
                end
            end
        end
    end

    return dmin, is, js

end

