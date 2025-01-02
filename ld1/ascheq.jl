#
# Copyright (C) 2004 PWSCF group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .

# Numerical integration of the radial schroedinger equation for
# bound states in a local potential.
# thresh determines the absolute accuracy for the eigenvalue
function ascheq!(nn, l, e, grid, vpot, Zval, thresh0, y, nstop)

    # nstop is a somewhat a return code, probably can be removed.
    # Its value depends on whether the calculation is successfull or not.

    Nrmesh = grid.Nrmesh
    println("input e = ", e)

    #
    #  set up constants and initialize
    #
    c = zeros(Float64, Nrmesh)
    f = zeros(Float64, Nrmesh)
    el = zeros(Float64, Nrmesh)

    thresh = thresh0
    
    if e < -5e2
        thresh = thresh0*10.0
    end
    
    ddx12 = grid.dx^2 / 12.0
    l1 = l + 1
    sqlhf = 0.5*(l + 0.5)^2 # Ha unit
    ndcr = nn - l - 1
    
    # set initial lower and upper bounds to the eigenvalue
    eup = vpot[Nrmesh] + sqlhf/grid.r2[Nrmesh]
    elw = eup
    for i in 1:Nrmesh
       elw = min( elw, vpot[i] + sqlhf/grid.r2[i] )
    end
    
    nstop = 200
    
    if eup == elw
        println("Should go to 900")
        #go to 900
    end
    
    if e > eup
        e = 0.9*eup + 0.1*elw
    end
    
    if e < elw
        e = 0.9*elw + 0.1*eup
    end

    println()
    println("Initial lower and upper bounds to the eigenvalue (in Ha)")
    println()
    @printf("e   = %18.10f\n", e)
    @printf("eup = %18.10f\n", eup)
    @printf("elw = %18.10f\n", elw)
    
    #
    #  series developement of the potential near the origin
    #
    for i in 1:4
        y[i] = vpot[i] - Zval/grid.r[i]
    end
    b = zeros(Float64,4) # XXX: b is originally b(0:3)
    radial_grid_series!( y, grid.r, grid.r2, b )

    println()
    println("Near origin: r, vpot, yi, b (in Ha)")
    println()
    for i in 1:4
        @printf("%18.10e %18.10e %18.10e %18.10e\n", grid.r[i], vpot[i], y[i], b[i])
    end

    #exit()

    NmaxIter = 50

    # needed outside the loop
    nstart = 0
    f2 = 0.0

    println()

    for iterSch in 1:NmaxIter
    
        println("===============================")
        println("iterSch = ", iterSch)
        println("===============================")

        nstop = 300
    
        #
        # set up the f-function and determine the position of its last
        # change of sign
        #
        # f < 0 (approximatively) means classically allowed   region
        # f > 0         "           "        "      forbidden   "
        #
        ik = 0
        # Using Numerov algorithm
        f[1] = ddx12*( grid.r2[1]*( vpot[1] - e ) + sqlhf)
        for i in 2:Nrmesh
            f[i] = ddx12*( grid.r2[i]*( vpot[i] - e ) + sqlhf)
            if f[i] != abs(f[i])*sign(f[i-1])
                ik = i
            end
        end
        println("f = ", f[1:4])
        nstop = 302
    
        # XXX: What's this?
        if ik >= (Nrmesh - 2)
            println("Line 106: Should go to 900")
            break
            #go to 900
        end

        for i in 1:Nrmesh
            f[i] = 1.0 - 2*f[i] # convert to Ry?
        end
        for i in 1:4
            @printf("f = %18.10f\n", f[i])
        end
        @printf("ddx12 = %18.10f\n", ddx12)
        fill!(y, 0.0)
  
        # determination of the wave-function in the first two points by
        # series developement
        xl1 = l + 1.0
        x4l6 = 4.0*l + 6.0
        b0e = b[1] - e # Ha unit
        c1 = Zval/xl1  # in Ha?
        c2 = (c1*Zval + b0e)/x4l6 # Ha
        println("e = ", e)
        start_scheq!( l, e, b, grid, Zval, y )
        @printf("After start_scheq! ")
        @printf("y[1] = %18.10f\n", y[1])
        @printf("y[2] = %18.10f\n", y[2])
        #if iter == 1
        #    println("exit ffr 138")
        #    exit()
        #end

        # start outward integration and count number of crossings
        ncross = 0
        ymx = 0.0
        println("\nStarting outward integration")
        println("ik = ", ik)
        println("fn = ", f[1:4])
        for n in 2:(ik-1)
            # Numerov algorithm
            y[n+1] = ( ( 12.0 - 10.0*f[n] )*y[n] - f[n-1]*y[n-1] ) / f[n+1]
            # Check for crossing here
            if y[n] != abs(y[n])*sign(y[n+1])
                ncross = ncross + 1
            end
            ymx = max( ymx, abs(y[n+1]) )
        end
        @printf("ymx = %18.10f\n", ymx)
        println("ncross = ", ncross)
        #if iter == 1
        #    println("exit ffr 159")
        #    exit()
        #end


        #  matching radius has been reached going out. if ncross is not
        #  equal to ndcr, modify the trial eigenvalue.
        if ndcr < ncross
            # too many crossings. e is an upper bound to the true eigenvalue.
            # increase abs(e)
            eup = e
            rap = ( Float64(ncross + l1)/nn )^2
            e = (e - vpot[Nrmesh] )*rap + vpot[Nrmesh]
            if e < elw
                e = 0.9*elw + 0.1*eup
            end
            # go to 300 # skip iter: continue
            continue
        elseif ndcr > ncross
            #  too few crossings. e is a lower bound to the true eigen-
            #  value. decrease abs(e)
            #
            elw = e
            rap = ( Float64(ncross+l1)/nn )^2
            e = ( e - vpot[Nrmesh] )*rap + vpot[Nrmesh]
            if e > eup
                e = 0.9*eup + 0.1*elw
            end
            #go to 300
            continue
        end
    
        println("ndcr = ", ndcr)
        #if iter == 1
        #    println("exit ffr 201")
        #    exit()
        #end


        # prepare inward integration
        # charlotte froese can j phys 41,1895(1963)
        #
        # start at  min(rmax, 10*rmatch)
        #
        nstart = Nrmesh
        ns = 10
        rstart = ns*grid.r[ik]
        if rstart < grid.r[Nrmesh]
            for i in ik:Nrmesh
                nstart = i
                if grid.r[i] >= rstart
                    #go to 403
                    break
                end
            end
            # 403 
            nstart = round(Int64, floor(nstart/2))
            nstart = 2*nstart + 1
        end
        println("nstart = ", nstart)
        #if iter == 1
        #    println("exit ffr 231")
        #    exit()
        #end


        # set up a, l, and c vectors
        n = ik + 1
        el[n] = 10.0*f[n] - 12.0
        c[n] = -f[ik]*y[ik]
        n2 = ik + 2
        for n in n2:nstart
           di = 10.0*f[n] - 12.0
           el[n] = di - f[n]*f[n-1]/el[n-1]
           c[n] = -c[n-1]*f[n-1]/el[n-1]
        end
        println("el[nstart] = ", el[nstart])
        println("c[nstart] = ", c[nstart])
        #if iter == 1
        #    println("exit ffr 247")
        #    exit()
        #end

        # Start inward integration by the Froese's tail procedure
        expn = exp( -sqrt( 12.0*abs(1.0 - f[nstart-1]) ) )
        y[nstart-1] = c[nstart-1]/( el[nstart-1] + f[nstart]*expn )
        y[nstart] = expn*y[nstart-1]
        for n in range(nstart-2, stop=ik+1, step=-1) #nstart-2,ik+1,-1
            y[n] = ( c[n] - f[n+1]*y[n+1])/el[n]
        end
        @printf("y = %18.10f\n", y[ik+1])
        #if iter == 1
        #    println("exit ffr 260")
        #    exit()
        #end

        # if necessary, improve the trial eigenvalue by the cooley's
        # procedure. jw cooley math of comp 15,363(1961)
        fe = ( 12.0 - 10.0*f[ik] )*y[ik] - f[ik-1]*y[ik-1] - f[ik+1]*y[ik+1]
        @printf("fe = %18.10f\n", fe)
        #if iter == 1
        #    println("exit ffr 270")
        #    exit()
        #end


        #  calculate the normalization
        if ymx >= 1.0e10
           for i in 1:Nrmesh
              y[i] = y[i]/ymx
           end
        end
        @printf("ymx = %18.10f\n", ymx)        
        @printf("y[1] = %18.10f\n", y[1])
        @printf("y[Nrmesh] = %18.10f\n", y[Nrmesh])
        #if iter == 1
        #    println("exit ffr 283")
        #    exit()
        #end
    
        a0 = 1.0/(2*l + 3)
        a1 = c1/(l + 2)
        a2 = (c1*c1 + c2 + c2)/(2*l + 5)
        sum0 = (a0 + grid.r[1]*( a1 + grid.r[1]*a2 ) )*grid.r[1]^(2*l + 3)
        nst2 = nstart - 2
        f2 = grid.r2[1]*y[1]*y[1]
        ss = grid.r[1]*f2/(2*l1 + 1)
        for n in range(1, stop=nst2, step=2)
            f0 = f2
            f1 = grid.r2[n+1]*y[n+1]*y[n+1]
            f2 = grid.r2[n+2]*y[n+2]*y[n+2]
            ss = ss + f0 + f2 + 4.0*f1
        end
        @printf("ss = %18.10f\n", ss)
        #if iter == 1
        #    println("exit ffr 310")
        #    exit()
        #end

        ss = sum0 + grid.dx*ss/3.0
        dfe = -y[ik]*f[ik]/grid.dx/ss
        de = -fe*dfe*0.5 # Ha unit?
        eeps = abs(de/e)
        
        @printf("ss = %18.10f\n", ss)
        @printf("de = %18.10f\n", de)

        #if iter == 4
        #    println("exit ffr 321")
        #    exit()
        #end


        println("iterSch = ", iterSch, " e = ", e,  " de = ", de)
        if abs(de) < 2*thresh
            println("GOTO 600 here ....")
            break
            #go to 600
            #continue
        end
    
        if eeps > 0.25
            println("Updating de in 341")
            de = 0.25*de/eeps
        end

        if de > 0.0
            println("Updating elw in 346")
            elw = e
        end
    
        if de < 0.0
            println("Updating eup")
            eup = e
        end
  
        e = e + de
        println("de = ", de)
        println("e = ", e)
  
        if e > eup
            println("Updating e 358")
            e = 0.9*eup + 0.1*elw
        end
    
        if e < elw
            println("Updating e 363")
            e = 0.9*elw + 0.1*eup
        end
    
        #if(iter .lt. maxter) go to 300
  
        nstop = 50
    
        println("New data:")
        @printf("e = %18.10f\n", e)
        @printf("eup = %18.10f\n", eup)
        @printf("elw = %18.10f\n", elw)

    end

    println("e = ", e)

    #println("exit ffr 340")
    #exit()

    #600 continue
    # LOOP


    # normalize the eigenfunction and exit
    for n in nstart:(Nrmesh-1)
        
        y[n+1] = 0.0
        
        if y[n] == 0.0 # FIXME
            continue
            #go to 601
        end
        
        yln = log(abs(y[n]))
        xp = -sqrt(12.0*abs(1.0 - f[n]))
        expn = yln + xp
        if expn < -80.0
            continue
            #go to 601
        end
        
        y[n+1] = abs(exp(expn))*sign(y[n])
        # 601 continue
    end


    sum1 = 0.0
    for n in range(nstart, stop=Nrmesh-2, step=2)
       f0 = f2
       f1 = grid.r2[n+1]*y[n+1]*y[n+1]
       f2 = grid.r2[n+2]*y[n+2]*y[n+2]
       sum1 = sum1 + f0 + f2 +4.0*f1
    end
    ss = ss + grid.dx*sum1/3.0
    ss = sqrt(ss)
    for n in 1:Nrmesh
        y[n] = grid.sqrtr[n]*y[n]/ss
    end
    
    if nstop < 100
        println("Should go to 900")
        return e, nstop
        #go to 900
    end
    
    nstop = 0
    
    return e, nstop

end
