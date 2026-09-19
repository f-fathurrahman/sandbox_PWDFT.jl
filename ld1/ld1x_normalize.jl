
function l1dx_normalize!(ld1x_input, grid, idx_rcut, qq, beta_prj, phi, ℓ)
    Nrmesh = grid.Nrmesh
    Nbeta = ld1x_input.Nbeta
    lls = ld1x_input.lls
    gi = zeros(Float64, Nrmesh)
    # if US pseudopotential compute the augmentation part
    nst = (ℓ + 1)*2
    work = zeros(Float64, Nbeta)
    for ibeta in 1:Nbeta
        if ℓ == lls[ibeta] # .and. abs(j-jjs(n1)) < 1.e-7_dp ) then
            idx_r = idx_rcut[ibeta]
            for ir in 1:idx_r
                gi[ir] = beta_prj[ir,ibeta]*phi[ir]
            end
            work[ibeta] = integ_0_inf_dr(gi, grid, idx_r, nst)
        else
            work[ibeta] = 0.0
        end
    end
    #
    for ir in 1:Nrmesh
        gi[ir] = phi[ir]*phi[ir]
    end
    work1 = integ_0_inf_dr(gi, grid, Nrmesh, nst)
    #
    # and adding to the charge density
    for ibeta in 1:Nbeta, jbeta in 1:Nbeta
        work1 += qq[ibeta,jbeta] * work[ibeta] * work[jbeta]  
    end

    if abs(work1) < 1e-10
        println("Zero norm: self consistency problem; state = $ns")
        work1 = 1.0
    elseif work1 <= -1e-10
        error("Negative norm")   
    end
    work1 = sqrt(work1)
    for ir in 1:Nrmesh
        phi[ir] = phi[ir]/work1
    end
    return
end
