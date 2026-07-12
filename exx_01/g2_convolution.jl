function g2_convolution!( exx, LatVecs, gvec_exx, xk, xkq, fac )

    Ng = gvec_exx.Ng
    G = gvec_exx.G
    SMALL = 1e-6
    SMALL_QDIV = 1e-8
    grid_factor = 1.0 #XXX should be from exx

    gau_scrlen = exx.gau_scrlen
    x_gamma_extrapolation = exx.x_gamma_extrapolation
    exxdiv = exx.exxdiv
    yukawa = exx.yukawa
    Nq1 = exx.Nq1
    Nq2 = exx.Nq2
    Nq3 = exx.Nq3

    q = zeros(Float64, 3)
    grid_factor_track = zeros(Float64, Ng)
    qq_track = zeros(Float64, Ng)
    odg = zeros(Bool, 3)

    # First the types of Coulomb potential that need q(3) and an external call
    #=
    if use_coulomb_vcut_ws
        for ig in 1:Ng
            q[:] = xk[:] - xkq[:] + G[:,ig]
            fac[ig] = vcut_get(vcut,q)
        end
        return
    end
    if use_coulomb_vcut_spheric
        for ig in 1:Ng
            q[:]= xk[:] - xkq[:] + G[:,ig]
            fac(ig) = vcut_spheric_get(vcut,q)
        end
        return
    end
    =#

    #
    # Now the Coulomb potential that are computed on the fly
    nqhalf = 0.5*[Nq1, Nq2, Nq3]
    #
    # Set the grid_factor_track and qq_track
    #
    if x_gamma_extrapolation
        for ig in 1:Ng
            q[:] = xk[:] - xkq[:] + G[:,ig]
            qq_track[ig] = sum(q.^2)
            #
            x = ( q[1]*LatVecs[1,1] + q(2)*LatVecs[2,1] + q(3)*LatVecs[3,1] ) * nqhalf[1]
            odg[1] = abs(x - round(Int64,x)) < SMALL
            #
            x = ( q[1]*LatVecs[1,2] + q[2]*LatVecs[2,2] + q[3]*LatVecs[3,2] ) * nqhalf[2]
            odg[2] = abs(x - round(Int64,x)) < SMALL
            #
            x = ( q[1]*LatVecs[1,3] + q[2]*LatVecs[2,3] + q[3]*LatVecs[3,3] ) * nqhalf[3]
            odg[3] = abs(x - round(Int64,x)) < SMALL
            #
            if all(odg)
                grid_factor_track[ig] = 0.0 # on double grid
            else
                grid_factor_track[ig] = grid_factor # not on double grid
            end
        end
    else
        # No gamma extrapolation
        for ig in 1:Ng
            q[:] = xk[:] - xkq[:] + G[:,ig]
            qq_track[ig] = sum(q.^2)
        end
        fill!(grid_factor_track, 1.0)
    end
    #
    # The big loop
    for ig in 1:Ng
        qq = qq_track[ig]
        #
        if gau_scrlen > 0
            fac[ig] = ( (pi/gau_scrlen)^1.5 )*exp(-qq/4/gau_scrlen) * grid_factor_track[ig]
        #
        elseif qq > SMALL_QDIV
            if erfc_scrlen > 0
                fac[ig] = 4*pi/qq*(1.0 - exp(-qq/4/erfc_scrlen^2)) * grid_factor_track[ig]
            elseif erf_scrlen > 0
                fac[ig] = 4*pi/qq*( exp(-qq/4/erf_scrlen^2) ) * grid_factor_track[ig]
            else
                fac[ig] = 4*pi/( qq + yukawa ) * grid_factor_track[ig]
            end
        #
        else
            #
            fac[ig] = -exxdiv # or rather something else (see F.Gygi)
            if yukawa > 0.0 && !x_gamma_extrapolation
                fac[ig] += 4*pi/( qq + yukawa ) # add additional factor
            end
            #      
            if erfc_scrlen > 0.0 && !x_gamma_extrapolation
                fac[ig] += pi/(erfc_scrlen^2)
            end
        end
    end
    return
end