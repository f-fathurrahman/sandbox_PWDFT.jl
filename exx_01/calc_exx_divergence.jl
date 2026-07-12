function calc_exx_divergence(pw, exx)
    
    if !use_regularization
        return 0.0
    end

    ecutwfc = pw.ecutwfc
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs

    gcutw = sqrt(2*ecutwfc)
    α = 10.0/gcutw

    Nq1 = exx.Nq1
    Nq2 = exx.Nq2
    Nq3 = exx.Nq3

    SMALL = 1e-6

    dq1 = 1.0/Nq1
    dq2 = 1.0/Nq2 
    dq3 = 1.0/Nq3 
    res = 0.0
    odg = zeros(Bool, 3)
    xq = zeros(Float64, 3)
    q = zeros(Float64, 3)
    on_double_grid = false
    for iq1 in 1:Nq1, iq2 in 1:Nq2, iq3 in 1:Nq3
        xq[:] = RecVecs[:,1] * (iq1-1) * dq1 +
                RecVecs[:,2] * (iq2-1) * dq2 +
                RecVecs[:,3] * (iq3-1) * dq3
        for ig in 1:Ng
            q[1] = xq[1] + G[1,ig]
            q[2] = xq[2] + G[2,ig]
            q[3] = xq[3] + G[3,ig]
            qq = ( q[1]^2 + q[2]^2 + q[3]^2 )
            if x_gamma_extrapolation
                x = ( q[1]*LatVecs[1,1] + q(2)*LatVecs[2,1] + q(3)*LatVecs[3,1] ) * nqhalf[1]
                odg[1] = abs(x - round(Int64,x)) < SMALL
                #
                x = ( q[1]*LatVecs[1,2] + q[2]*LatVecs[2,2] + q[3]*LatVecs[3,2] ) * nqhalf[2]
                odg[2] = abs(x - round(Int64,x)) < SMALL
                #
                x = ( q[1]*LatVecs[1,3] + q[2]*LatVecs[2,3] + q[3]*LatVecs[3,3] ) * nqhalf[3]
                odg[3] = abs(x - round(Int64,x)) < SMALL
                #
                on_double_grid = all(odg)
            end
            if !on_double_grid
                if qq > 1e-8
                    if erfc_scrlen > 0
                        res += exp(-α*qq) / qq * (1.0 - exp(-pi^2*qq/erfc_scrlen^2)) * grid_factor
                    elseif erf_scrlen >0
                        res += exp(-α*qq) / qq * (exp(-pi^2*qq/erf_scrlen^2)) * grid_factor
                    else
                        res += exp(-α*qq) / (qq + yukawa/(4*pi^2)) * grid_factor
                    end
                end
            end
        end
    end
    
    if !x_gamma_extrapolation
        if yukawa > 0.0
            res += 4*pi^2/yukawa
        elseif erfc_scrlen > 0.0
            res += pi^2/erfc_scrlen^2
        else
            res -= α
        end
    end
    res = res*pi/nqs
    α = α / (4*pi^2)
    nqq = 100000
    dq = 5.0 / SQRT(α) / nqq
    aa = 0.0
    for iq in 0:nqq
        q_ = dq * (iq + 0.5d0)
        qq = q_ * q_
        if erfc_scrlen > 0
            aa -= exp(-α*qq) * exp(-qq/4/erfc_scrlen^2)*dq
        elseif erf_scrlen > 0
            aa = 0.0
        else
            aa -= -exp(-α*qq) * yukawa / (qq + yukawa)*dq
        end
    end
    aa *= 8/(4*pi)
    aa += 1.0/sqrt(α*π)
    if erf_scrlen > 0
        aa = 1.0/sqrt( (α + 1.0/4.0/erf_scrlen^2) * pi )
    end
    res -= CellVolume*aa
    return res*nqs
end