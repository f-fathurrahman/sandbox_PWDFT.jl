function seriesbes!(f, r, r2, npt, xc)
    npt2 = round(Int64, npt/2+1)
    xc[3] = ( (f[1] -f[npt2])/(r[1]-r[npt2]) -(f[npt]-f[npt2])/(r[npt]-r[npt2]) )/ ( r[1]-r[npt] )
    xc[1] = f[1]
    xc[2] = ( f[npt]-f[npt2] ) / (r[npt]-r[npt2]) -xc[3]*(r[npt]+ r[npt2] )
    xc[4] = 0.0
    return
end


function _solve_tridiag!(a, b, c, r, u, n)
    gam = zeros(Float64, n)

    if abs(b[1]) < 1e-10
        error("b[1] is too small")
    end
    bet = b[1]
    u[1] = r[1]/bet
    for j in 2:n
        gam[j] = c[j-1]/bet
        bet = b[j] - a[j]*gam[j]
        if abs(bet) < 1e-10
            error("bet is too small")
        end
        u[j] = ( r[j] - a[j]*u[j-1] )/bet
    end
    for j in range(n - 1, stop = 1, step = -1)
        u[j] = u[j] - gam[j+1]*u[j+1]
    end
    return
end

function compute_chi!(grid, V_Ps_loc, ℓ, idx_r, phi_in, chi_out, xc, E_in, lbes4)
    #= 
    This routine computes the chi functions:
        |chi> = (\epsilon - T - V_{loc}) |psi>
    =#

    Nrmesh = grid.Nrmesh
    #println("compute_chi: idx_r = $idx_r")

    j1 = zeros(Float64, Nrmesh)
    aux = zeros(Float64, Nrmesh)
    gi = zeros(Float64, Nrmesh)
    b = zeros(Float64, 4)
    c = zeros(Float64, 4)
    arow = zeros(Float64, Nrmesh)
    brow = zeros(Float64, Nrmesh)
    crow = zeros(Float64, Nrmesh)
    drow = zeros(Float64, Nrmesh)
    
    # RRKJ: first expand in a taylor series the phis function
    # Since we know that the phis functions are a sum of Bessel 
    # functions with coefficients xc, we compute analytically
    # the asymptotic expansion
    ddx12 = grid.dx^2 / 12.0
    x4l6 = 4*ℓ + 6
    x6l12 = 6*ℓ + 12

    for i in 1:6
        j1[i] = phi_in[i]/grid.r[i]^(ℓ + 1)
    end
    seriesbes!(j1, grid.r, grid.r2, 6, c)
    
    if ℓ == 0
        if lbes4 || rho0 == 0.0
            c[1] = xc[1] + xc[2] + xc[3]
            c[2] = 0.0
            c[3] = -xc[1]*(xc[4]^2/6.0) - xc[2]*(xc[5]^2/6.0) - xc[3]*(xc[6]^2/6.0)
            c[4] = 0.0
        else
            c[1] = xc[1]+xc[2]+xc[3]+xc[4]
            c[2] = 0.0
            c[3] = -xc[1]*(xc[5]^2/6.0) -xc[2]*(xc[6]^2/6.0) - xc[3]*(xc[7]^2/6.0) - xc[4]*(xc[8]^2/6.0)
            c[4] = 0.0
        end
    #
    elseif ℓ == 3
        c[1] = xc[1]*(48.0*xc[4]^3/5040.0) + xc[2]*(48.0*xc[5]^3/5040.0) + xc[3]*(48.0*xc[6]^3/5040.0)
        c[2] = 0.0
        c[3] = -xc[1]*(192.0*xc[4]^5/362880.0) - xc[2]*(192.0*xc[5]^5/362880.0) - xc[3]*(192.0*xc[6]^5/362880.0)
        c[4] = 0.0
    #
    elseif ℓ == 2
        c[1] = xc[1]*(xc[4]^2/15.0) + xc[2]*(xc[5]^2/15.0) + xc[3]*(xc[6]^2/15.0)
        c[2] = 0.0
        c[3] = -xc[1]*(xc[4]^4/210.0) - xc[2]*(xc[5]^4/210.0) - xc[3]*(xc[6]^4/210.0)
        c[4] = 0.0
    #
    elseif ℓ == 1
        c[1] = xc[1]*(xc[4]/3.0) + xc[2]*(xc[5]/3.0) + xc[3]*(xc[6]/3.0)
        c[2] = 0.0
        c[3] = -xc[1]*(xc[4]^3/30.0) - xc[2]*(xc[5]^3/30.0) - xc[3]*(xc[6]^3/30.0)
        c[4] = 0.0
    else
        error("ℓ=$ℓ is not programmed")
    end
  
    # and the potential
    for i in 1:4
        j1[i] = V_Ps_loc[i]
    end
  
    if abs( j1[1] - j1[4]) > 1e-12
        #call series(j1, grid%r, grid%r2, b)
        radial_grid_series!( j1, grid.r, grid.r2, b )
    else
        b = 0.0
        b[1] = j1[1]
    end
    #
    # and compute the taylor expansion of the chis
    b0e = (b[1] - E_in)*2.0 # this is in Hartree
    g0 = x4l6*c[3] - b0e*c[1]
    g1 = x6l12*c[4] - c[1]*b[2]
    g2 = -(b0e*c[3] + b[3]*c[1])
    ir_start = 5
    for ir in 1:ir_start-1
        chi_out[ir] = (g0 + grid.r[ir]*(g1 + g2*grid.r[ir]))*grid.r[ir]^(ℓ+3)/grid.sqrtr[ir]
    end
    #
    for ir in 1:Nrmesh
        aux[ir] = (g0 + grid.r[ir]*(g1 + g2*grid.r[ir]))
    end
    #
    # set up the equation
    for ir in 1:Nrmesh
        gi[ir] = phi_in[ir]/grid.sqrtr[ir]
    end
    for ir in 1:Nrmesh
        j1[ir] = 2*grid.r2[ir]*(V_Ps_loc[ir] - E_in) + (ℓ + 0.5)^2 # this is in Hartree
        j1[ir] = 1.0 - ddx12*j1[ir]
    end
    #
    for ir in ir_start:Nrmesh-3
        drow[ir] = gi[ir+1]*j1[ir+1] + gi[ir]*(-12.0 + 10.0*j1[ir]) + gi[ir-1]*j1[ir-1]
        brow[ir] = 10*ddx12
        crow[ir] = ddx12
        arow[ir] = ddx12
    end
    drow[ir_start] = drow[ir_start] - ddx12*chi_out[ir_start-1]
    chi_out[Nrmesh-2] = 0.0
    chi_out[Nrmesh-1] = 0.0
    chi_out[Nrmesh] = 0.0
    #
    # and solve it
    Nvars = Nrmesh - 3 - ir_start
    idx = ir_start:(ir_start+Nvars-1)
    @views _solve_tridiag!( arow[idx], brow[idx], crow[idx], drow[idx], chi_out[idx], Nvars)
    #
    # put the correct normalization and r dependence
    #  
    for ir in 1:Nrmesh
        chi_out[ir] = chi_out[ir] * grid.sqrtr[ir]/grid.r2[ir]
    end
    #
    # smooth close to the origin with asymptotic expansion
    for ir in ir_start:Nrmesh
        if abs(chi_out[ir]/grid.r[ir]^(ℓ+1) - aux[ir] ) < 1e-3*abs(aux[ir])
            break #goto 100 ! break
        end
        chi_out[ir] = aux[ir]*grid.r[ir]^(ℓ+1)
    end

    # LABEL100
    #if n == grid%mesh+1 .or. grid%r(min(n,grid%mesh)) > 0.05)then
    #write(stdout,*) lam,n,grid%mesh,grid%r(min(n,grid%mesh))
    #call errore('compute_chi','n is too large',1)
    #endif
    #
    # When the input wavefunction is a diverging scattering state 
    # this routine might become numerically unstable at large r. 
    # Here chi_out should be 0.0 but it is not due to this instability. 
    # We now clean chi_out when phi_in > 20.0.
    # Clean also after 7.0 a.u..
    #     
    r_clean = 100.0
    for ir in 1:Nrmesh
        if abs(phi_in[ir]) > 20.0
            r_clean = grid.r[ir]
            break
        end
    end
    println("r_clean = ", r_clean)
    r_clean = min(r_clean, 7.0)
    
    if r_clean < grid.r[idx_r]
        error("phi_in too large before r_c")
    end
    
    for ir in range(Nrmesh, stop = 1, step = -1)
        if grid.r[ir] < r_clean
            break
        end
        chi_out[ir] = 0.0
    end

    # LABEL200
    # check that the chi are zero beyond ikk
    nst = 0
    fill!(gi, 0.0)
    for ir in (idx_r+1):Nrmesh
        gi[ir] = chi_out[ir]^2
    end
    for ir in min(idx_r+20, grid.Nrmesh):Nrmesh
        chi_out[ir] = 0.0
    end
    integ = integ_0_inf_dr(gi, grid, Nrmesh, nst)
    verbose = false
    if integ > 2.e-6
        println("ℓ=$ℓ integ=$integ rcut=$(grid.r[idx_r])")
        if verbose
            for ir in idx_r:Nrmesh
                println("$(grid.r[ir]) $(gi[ir])")
            end
        end
        @warn("chi too large beyond r_c")
    end
    return
end