function calc_charge_ps!(rho_i, phi_i, nwf_i, ll_i, jj_i, oc_i, iswf_i)

    # HARCODED
    which_augfun = :PSQ

    work = zeros(Float64, nwf_i)
    gi = zeros(Float64, Nrmesh)

    fill!(rho_i, 0.0)
    #
    # compute the square modulus of the eigenfunctions
    for iwf in 1:nwf_i
        if oc_i[iwf] > 0.0
            ispin = iswf_i[iwf]
            for ir in 1:Nrmesh
                rho_i[ir,ispin] += oc_i[iwf] * phi_i[ir,iwf]^2
            end
        end
    end
    #
    # if US pseudopotential compute the augmentation part
    #
    #if( pseudotype == 3 ) then
    # XXXWe assume that all pseudotype is 3 ?
    for iwf in 1:nwf_i
        # skip is occupation is zero or negative
        if oc_i[iwf] <= 0.0
            continue
        end
        ispin = iswf_i[iwf]
        for ibeta in 1:Nbeta
            if ll_i[iwf] == lls[ibeta] #.and. abs(jj_i(ns)-jjs(n1)) < 1.e-7_dp) then
                nst = (ll_i[iwf] + 1)*2
                idx_r = idx_rbeta[ibeta]
                for ir in 1:idx_r
                    gi(n) = betas(n,n1)*phi_i(n,ns)
                end
                work[ibeta] = int_0_inf_dr(gi, grid, idx_r, nst)
            else
                work[ibeta] = 0.0
            end
        end
        #
        # and adding to the charge density
        for ibeta in 1:Nbeta, jbeta = 1,nbeta
            if which_augfun == :PSQ
                for ir in 1:Nrmesh
                    rho_ip[ir,ispin] += qvanl[ir,ibeta,jbeta,0] * oc_i[iwf] * work[ibeta] * work[jbeta]
                end
            else
                for ir in 1:Nrmesh
                    rho_i[ir,ispin] += qvan[ir,ibeta,jbeta] * oc_i[iwf] * work[ibeta] * work[jbeta]
                end
            end
        end
    end
  
    #
    # Check for negative charge
    for ispin in 1:Nspin
        for ir in 2:Nrmesh # ffr: why start from 2 ?
            if rho_i[ir,is] < -1.d-12
                error("negative rho found")
            end
        end
    end

    return
end