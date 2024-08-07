function solve_atom!(
    ptnucl::Bool, zn, nst, n, l, k, occ,
    xc_calc, xcgrad, r, evals, rho, vr, rwf;
    sol=137.035999084, maxscl=200
)

    # potential convergence tolerance
    SMALL = 1.0e-6

    @assert nst > 0

    nr = size(r, 1)

    # allocate local arrays
    vn = zeros(Float64,nr)
    vh = zeros(Float64,nr)
    vxc = zeros(Float64,nr)
    vrp = zeros(Float64,nr)
    
    ri = zeros(Float64,nr)
    wpr = zeros(Float64,4,nr)
    fr1 = zeros(Float64,nr)
    fr2 = zeros(Float64,nr)
    gr1 = zeros(Float64,nr)
    gr2 = zeros(Float64,nr)
    
    #if (xcgrad.eq.1) then
    #  allocate(grho(nr),g2rho(nr),g3rho(nr))
    #end if

    # find total electronic charge
    ze = 0.0
    for ist in 1:nst
        ze = ze + occ[ist]
    end

    # set up nuclear potential
    # FIXME: use vcnl from atsp_vars?
    potnucl!(ptnucl, r, zn, vn)

    for ir in 1:nr
        ri[ir] = 1.0/r[ir]
        # initialize the Kohn-Sham potential to the nuclear potential
        vr[ir] = vn[ir]
    end

    # determine the weights for radial integration
    wsplintp!(nr, r, wpr)
    
    # initialize mixing parameter
    betamix = 0.5
    
    # initialize eigenvalues to relativistic values (minus the rest mass energy)
    for ist in 1:nst
        t1 = sqrt( k[ist]^2 - (zn/sol)^2 )
        t1 = ( n[ist] - abs(k[ist]) + t1)^2
        t1 = 1.0 + (zn/sol)^2/t1
        evals[ist] = sol^2/sqrt(t1) - sol^2
    end
    
    dvp = 0.0
    # start of self-consistent loop
    is_converged = false
    for iscl in 1:maxscl

        # solve the Dirac equation for each state
        for ist in 1:nst
            @views evals[ist] = rdirac!( n[ist], l[ist], k[ist], r, vr, evals[ist],
                rwf[:,1,ist], rwf[:,2,ist]; sol=sol )
        end
    
        # compute the charge density
        for ir in 1:nr
            ss = 0.0
            for ist in 1:nst
              ss = ss + occ[ist]*( rwf[ir,1,ist]^2 + rwf[ir,2,ist]^2 )
            end
            fr1[ir] = ss
            fr2[ir] = ss*ri[ir]
            rho[ir] = ss*ri[ir]^2 / (4*pi)
        end
        splintwp!(nr, wpr, fr1, gr1)
        splintwp!(nr, wpr, fr2, gr2)
        
        # find the Hartree potential
        t1 = gr2[nr]
        for ir in 1:nr
           vh[ir] = gr1[ir]*ri[ir] + t1 - gr2[ir]
        end
    
        # normalize charge density and potential
        t1 = ze/gr1[nr]
        for ir in 1:nr
            rho[ir] = t1*rho[ir]
            vh[ir] = t1*vh[ir]
        end
    
        # compute the exchange-correlation energy and potential
        if xcgrad == 1
            # GGA functional
            # |grad rho|
            fderiv!(1, nr, r, rho, grho)
            # grad^2 rho
            fderiv!(2, nr, r, rho, g2rho)
            for ir in 1:nr
                g2rho[ir] = g2rho[ir] + 2.0*ri[ir]*grho[ir]
            end
            # approximate (grad rho).(grad |grad rho|)
            for ir in 1:nr
                g3rho[ir] = grho[ir]*g2rho[ir]
            end
            #xcifc(xctype,n=nr,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, ec=ec,vx=vx,vc=vc)
        else
            # LDA functional
            calc_Vxc_LDA!( xc_calc, rho, vxc )
        end
    
        # self-consistent potential
        for ir in 1:nr
            vr[ir] = vh[ir] + vxc[ir]
        end
    
        # determine change in potential
        ss = 0.0
        for ir in 1:nr
          ss = ss + ( vr[ir] - vrp[ir] )^2
        end
        dv = sqrt(ss)/nr
    
        if iscl > 2
            # reduce beta if change in potential is diverging
            if dv > dvp
                betamix = betamix*0.8
            end
            betamix = max(betamix,0.01)
        end
    
        dvp = dv
    
        for ir in 1:nr
            # mix old and new potentials
            vr[ir] = (1.0 - betamix)*vrp[ir] + betamix*vr[ir]
            vrp[ir] = vr[ir]
            # add nuclear potential
            vr[ir] = vr[ir] + vn[ir]
        end

        #@printf("iscl = %5d dv = %10.5e\n", iscl, dv)

        # check for convergence
        if ( (iscl > 2) && (dv < SMALL) )
            println(".... Converged ....")
            is_converged = true
            break
        end
    end

    if !is_converged
        println("Warning(atom): maximum iterations exceeded")
    end

    return

end 

