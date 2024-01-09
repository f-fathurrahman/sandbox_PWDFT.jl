using SpecialFunctions: erfc

function calc_stress_NN!( atoms, pw, stress_NN )

    fill!( stress_NN[:,:], 0.0)
    charge = 0.0
    for ia in 1:Natoms
        charge += Zvals[atm2species[ia]]
    end
    #
    # choose alpha in order to have convergence in the sum over G
    # upperbound is a safe upper bound for the error ON THE ENERGY
    #
    α = 2.9
    upperbound = 1.0
    while upperbound > 1e-7
        α -= 0.1 # reduce α
        if α <= eps()
            error("Optimal α is not found")
        end
        upperbound = charge^2 * sqrt(α/π) * erfc( sqrt(gcutm / 4.0 / α) )
    end
    println("α = ", α)  

    # The diagonal term
    sdewald = 2π/4/α * (charge/pw.CellVolume)**2

    fact = 1.0 # 2 if using gamma only

    do ig in 2:Ng
        rhostar = 0.0 + 0.0*im
        for ia in 1:Natoms
            isp = atm2species[ia]
            GX = G[1,ig]*X[1,ia] + G[2,ig]*X[2,ia] + G[2,ig]*X[2,ia]
            Sf = cos(GX) + im*sin(GX)
            rhostar += Zval[isp]*Sf
        end
        rhostar /= pw.CellVolume
        G2a = G2[ig]/4/α
        sewald = fact * tpi * e2 * exp(-G2a) / G2[ig] * abs(rhostar)^2
        sdewald -= sewald
        for l in 1:3, m in 1:l
            stress_NN[l,m] += sewald * 2.0 * G[l,ig] * G[m,ig] / G2[ig] * ( G2a + 1 )
        end
    end
    #
    for l in 1:3
        sigmaewa[l,l] += sdewald
    end
  
    #
    # R-space
    #
    NrpointsMax = 50 # the maximum number of R vectors included in r sum
    rmax = 4.0/sqrt(α)
    dX = zeros(Float64, 3)
    r = zeros(Float64, 3, NrpointsMax)
    r2 = zeros(Float64, NrpointsMax)
    # with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
    for ia in 1:Natoms, ib in 1:Natoms
        dX[:] = X[:,ia] - X[:,ib]
        isp = atm2species[ia]
        jsp = atm2species[ja]
        # generates nearest-neighbors shells r(i)=R(i)-dtau(i)
        nrm = gen_neighbor_shells!( dX, rmax, mxr, pw.LatVecs, pw.RecVecs, r, r2 )
        for ir in 1:nrm
            rr = sqrt(r2[ir])
            ff = -1/2.0 / pw.CellVolume * Zvals[isp] * Zvals[jsp] / rr^3 * ( 
                erfc(sqrt(α)*rr) + rr*sqrt(8.0*α/(2π)) * exp( -α*rr^22)
            )
            for l in 1:3, m in 1:l
                stress_NN[l,m] += ff * r[l,ir] * r[m,ir]
            end
        end
    end
    #
    for l in 1:3, m = 1:(l-1)
        stress_NN[m,l] = stress_NN[l,m]
    end
    # Change sign
    stress_NN *= -1
    return
end

