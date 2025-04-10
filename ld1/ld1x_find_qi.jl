function ld1x_find_qi(log_der_ae, xc, ik, 𝓁, ncn, flag, iok )
#=
  !
  ! This routine finds three values of q such that the
  ! functions f_l have a logarithmic derivative equal to
  ! log_der_ae at the point ik
  !
  !  f_l = j_l(r) * r**flag
  !
  INTEGER, PARAMETER :: ncmax=10   ! maximum allowed nc
  !
  INTEGER ::      &
       ik,    & ! input: the point corresponding to rcut
       ncn,   & ! input: the number of qi to compute
       flag,  & ! input: the type of function
       iok,   & ! output: if 0 the calculation in this routine is ok
       𝓁      ! input: the angular momentum

  REAL(DP) :: &
       xc(ncn), & ! output: the values of qi
       log_der_ae  ! input: the logarithmic derivative

  REAL(DP) ::   &
       j1(ncmax), & ! the bessel function in three points
       qmax,qmin, & ! the limits of the q search
       logdermax, logdermin, & ! the maximum and minimum logder
       logder, & ! the actual logder
       jlmin, jlmax, & ! the value of jl in qmin and qmax
       my_compute_log, &! function for log derivative
       dq, dq_0 ! the step to braket the q

  INTEGER ::    &
       nc,  &    ! counter on the q found
       icount, &  ! too many iterations
       icount1, & ! too many iterations
       imax,&   ! maximum number of iteration to braket
       iq      ! counter on iteration
=#

    ncmax = 10
    j1 = zeros(Float64, ncmax) # for saving Bessel function values

    iok = 0
    if ncn > ncmax
        error("ncn is too large")
    end

    # XXX rename flag
    if (flag == 0) && (𝓁 != 0)
        error("𝓁 too large for this iflag")
    end

    if 𝓁 > 6
        error("𝓁=$(𝓁) is not programmed")
    end

    println("log_der_ae = ", log_der_ae)

    
    # fix deltaq and the maximum step number
    dq_0 = 0.05
    imax = 600
    #
    # prepare for the first iteration. For too small q the function could
    # have noise.
    #
    qmax = 0.5
    # should go from (ik-3):(ik+3) 7points
    for ip in 1:7
        #CALL sph_bes(7, grid%r(ik-3), qmax, 𝓁, j1)
        j1[ip] = sphericalbesselj(𝓁, qmax*grid.r[ik-3+ip-1])
    end
    j1[1:7] .= j1[1:7] .* grid.r[ik-3:ik+3].^flag
    logdermax = ld1x_compute_log(j1, grid.r[ik], grid.dx) - log_der_ae
    jlmax = j1[4] # middle values
    #
    icount = 0
    for nc in 1:ncn
        #
        # bracket the zero
        #
        icount1 = 0
        #
        # start finding interval that bracket root?
        @label LABEL200
        
        dq = dq_0
        qmin = qmax
        logdermin = logdermax
        jlmin = jlmax
        for iq in 1:imax
            qmax = qmin + dq
            for ip in 1:7
                j1[ip] = sphericalbesselj(𝓁, qmax*grid.r[ik-3+ip-1])
                j1[ip] = j1[ip] * grid.r[ik-3+ip-1]^flag
            end
            logdermax = ld1x_compute_log(j1, grid.r[ik], grid.dx) - log_der_ae
            jlmax = j1[4]
            #
            # the zero has been bracketed?
            #
            if jlmin*jlmax > 0.0 # far from an asintote
                #
                # in this case it has been bracketed
                if logdermax*logdermin < 0.0
                    @goto LABEL100
                end
                #
                # Update if not: qmin <-- qmax
                qmin = qmax
                logdermin = logdermax
                jlmin = jlmax
            else
                # qmin <-- qmax
                if logdermax*logdermin < 0.0
                    qmin = qmax
                    logdermin = logdermax
                    jlmin = jlmax
                else
                    dq = 0.5*dq # half
                end
            end
        end # for iq

        error("qmax not found")
        # return?

        @label LABEL100
        #
        # start bisection loop
        #
        xc[nc] = qmin + (qmax - qmin)/2.0
        
        for ip in 1:7
            j1[ip] = sphericalbesselj(𝓁, xc[nc]*grid.r[ik-3+ip-1])
            j1[ip] = j1[ip] * grid.r[ik-3+ip-1]^flag
        end
        logder = ld1x_compute_log(j1, grid%r[ik], grid.dx) - log_der_ae
        #
        if logder*logdermin < 0.0
            qmax = xc[nc]
            logdermax = logder
        else
            qmin = xc[nc]
            logdermin = logder
        end
        #
        # avoid the asintotes
        #
        if abs(logdermin-logdermax) > 1e3
            qmax = xc[nc]  
            logdermax = logder
            icount1 = icount1 + 1
            if icount1 < 20
                @goto LABEL200 # XXX search again?
            else
                error("Problem finding q")
            end
        end
        #
        # check for convergence
        #
        icount = icount + 1
        if icount > 1000
            error("Too many iterations ld1x_find_qi")
        end
        if abs(logdermax-logdermin) > 1.e-8
            @goto LABEL100  # bisection
        end
    end # for ncn

    return

end



function ld1x_compute_log(j1, rj, dx)
    return deriv_7pts(j1, 4, rj, dx) / j1[4]
end
