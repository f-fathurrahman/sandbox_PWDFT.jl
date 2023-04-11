function PAW_rho_lm!( i, becsum, pfunc, rho_lm, aug )
    # !! Sum up pfuncs x occupation to build radial density's angular momentum components.
    # !
    # USE ions_base,         ONLY : nat 
    # USE lsda_mod,          ONLY : nspin
    # USE noncollin_module,  ONLY : nspin_mag
    # USE uspp_param,        ONLY : upf, nh, nhm
    # USE uspp,              ONLY : indv, ap, nhtolm,lpl,lpx
    # USE constants,         ONLY : eps12
    # USE atom,              ONLY : g => rgrid
    # !
    # TYPE(paw_info), INTENT(IN) :: i
    # !! atom's minimal info
    # REAL(DP), INTENT(IN)  :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
    # !! cross band occupation
    # REAL(DP), INTENT(IN)  :: pfunc(i%m,i%b,i%b)
    # !! psi_i * psi_j
    # REAL(DP), INTENT(OUT) :: rho_lm(i%m,i%l**2,nspin_mag)
    # !! AE charge density on rad. grid
    # REAL(DP), OPTIONAL,INTENT(IN) :: aug(i%m,(i%b*(i%b+1))/2,0:2*upf(i%t)%lmax)
    # !! augmentation functions (only for PS part)
    # !
    # ! ... local variables
    # !
    # REAL(DP) :: pref ! workspace (ap*becsum)
    # !
    # INTEGER :: ih, jh, &      ! counters for pfunc ih,jh = 1, nh (CRYSTAL index)
    #         nb, mb, &      ! counters for pfunc nb,mb = 1, nbeta (ATOMIC index)
    #         ijh, nmb,  &   ! composite "triangular" index for pfunc nmb = 1,nh*(nh+1)/2
    #         lm, lp, l, &   ! counters for angular momentum lm = l**2+m
    #         ispin          ! counter for spin (FIXME: may be unnecessary)
  
    # initialize density
    fill!(rho_lm, 0.0)
    
    for ispin in 1:Nspin
        ijh = 0 
        # loop on all pfunc for this kind of pseudo 
        for ih in 1:nh[isp], jh in ih:nh[isp]
            ijh = ijh + 1 
            nb = indv[ih,isp] 
            mb = indv[jh,isp] 
            nmb = (mb*(mb-1))/2 + nb  # mb has to be >= nb 
            if abs(becsum[ijh,ia,ispin]) < 1e-12
                continue
            end 
            # Loop over angular momentum angular_momentum: & 
            for lp in 1:lpx[nhtolm[jh,isp], nhtolm[ih,isp]] # lmaxq**2 
                # the lpl array contains the possible combination of LM,lm_j,lm_j that 
                # have non-zero a_{LM}^{(lm)_i(lm)_j} (it saves some loops) 
                lm = lpl(nhtolm(jh,i%t), nhtolm(ih,i%t), lp) 
                # 
                # becsum already contains a factor 2 for off-diagonal pfuncs 
                pref = becsum[ijh,ia,ispin] * ap[lm, nhtolm[ih,isp], nhtolm[jh,isp]] 
                #
                rho_lm[1:Nrmesh,lm,ispin] .+= pref * pfunc[1:Nrmesh,nb,mb] 
                if aug != nothing
                    # if I'm doing the pseudo part I have to add the augmentation charge 
                    l = floor(Int64, sqrt(lm-1) ) # l has to start from zero, lm = l**2 +m 
                    rho_lm[1:Nrmesh,lm,ispin] .+= pref * aug[1:Nrmesh,nmb,l] 
                end
            end # lp  
        end # ih, jh
    end
    return
end