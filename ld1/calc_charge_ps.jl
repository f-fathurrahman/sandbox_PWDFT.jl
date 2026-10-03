function calc_charge_ps!(
    grid, rho_i, phi_i, nwf_i, ll_i, oc_i,
    beta_prj, Nbeta, lls, idx_rbeta, qvan, qvanl;
    Nspin = 1
)

    # HARCODED
    which_augfun = :PSQ
    ispin = 1
    @assert Nspin == 1

    Nrmesh = grid.Nrmesh

    @info "Nbeta, nwf_i: " Nbeta nwf_i

    work = zeros(Float64, Nbeta)
    gi = zeros(Float64, Nrmesh)

    fill!(rho_i, 0.0)
    #
    # compute the square modulus of the eigenfunctions
    for iwf in 1:nwf_i
        if oc_i[iwf] > 0.0
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

        println("\niwf=$iwf oc_i=$(oc_i[iwf])")
        # skip is occupation is zero or negative
        if oc_i[iwf] <= 0.0
            println("This iwf=$iwf is skipped")
            continue
        end
        for ibeta in 1:Nbeta
            if ll_i[iwf] == lls[ibeta] #.and. abs(jj_i(ns)-jjs(n1)) < 1.e-7_dp) then
                nst = (ll_i[iwf] + 1)*2
                idx_r = idx_rbeta[ibeta]
                for ir in 1:idx_r
                    gi[ir] = beta_prj[ir,ibeta] * phi_i[ir,iwf]
                end
                work[ibeta] = integ_0_inf_dr(gi, grid, idx_r, nst)
            else
                work[ibeta] = 0.0
            end
            println("ibeta = $ibeta work = $(work[ibeta])")
        end
        #
        # and adding to the charge density
        for ibeta in 1:Nbeta, jbeta in 1:Nbeta
            #
            fij = abs(oc_i[iwf] * work[ibeta] * work[jbeta])
            #if fij > 0
            #    println("should contribute: ibeta=$ibeta jbeta=$jbeta fij=$fij")
            #end
            #
            if which_augfun == :PSQ
                for ir in 1:Nrmesh
                    rho_i[ir,ispin] += qvanl[ir,ibeta,jbeta,0] * fij
                end
            else
                for ir in 1:Nrmesh
                    rho_i[ir,ispin] += qvan[ir,ibeta,jbeta] * fij
                end
            end
        end
    end
  
    #
    # Check for negative charge
    for ispin in 1:Nspin
        for ir in 2:Nrmesh # ffr: why start from 2 ?
            if rho_i[ir,ispin] < -1e-12
                error("negative rho found")
            end
        end
    end

    return
end