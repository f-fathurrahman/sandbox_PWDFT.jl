# This is Julia version of hydrogen_radial.f90 program written by P. Giannozzi
# Added alternative way to match the outward and inward wave function

using Printf
import Plots, PlotThemes
Plots.theme(:dark)

function initialize_mesh(zmesh, xmin, dx, Nrmesh)
    r = zeros(Float64, Nrmesh)
    sqr = zeros(Float64, Nrmesh)
    r2 = zeros(Float64, Nrmesh)

    for i in 1:Nrmesh
        x = xmin + dx * (i-1)
        r[i] = exp(x) / zmesh
        sqr[i] = sqrt(r[i])
        r2[i] = r[i] * r[i]
    end

    @printf("\n radial grid information:\n")
    @printf(" dx =%10.6f, xmin =%10.6f, zmesh =%10.6f\n", dx, xmin, zmesh)
    @printf(" Nrmesh =%6d, r0=r[1] =%10.6f, r(Nrmesh) =%10.6f\n", Nrmesh, r[1], r[end])
    println()

    return r, sqr, r2
end

function initialize_V_atom(Zatom, r, Nrmesh)
    vpot = zeros(Float64, Nrmesh)
    for i in 1:Nrmesh
        vpot[i] = -Zatom/r[i]
    end
    return vpot
end

function solve_sch_rad(n, l, Nrmesh, dx, r, sqr, r2, vpot, Zatom)
    maxiter = 100
    TOL = 1.0e-10

    f = zeros(Float64, Nrmesh)
    y = zeros(Float64, Nrmesh)

    ddx12 = dx * dx / 12.0
    sqlhf = (l + 0.5)^2
    x2l2 = 2 * l + 2

    eup = vpot[end]
    elw = minimum( sqlhf ./ (2.0 .* r2) .+ vpot )
    if eup - elw < TOL
        println(elw, " ", eup)
        error("ERROR: solve_sch_rad: lower and upper bounds are equal")
    end

    E = 0.5 * (elw + eup)
    de = 0.0
    converged = false
    iter = 0

    nodes = n - l - 1
    ncross = -1
    icl = -1

    E_old = Inf

    for kkk in 1:maxiter
        iter = kkk

        # Set up f-function and find the last sign change (classical turning point)
        icl = -1
        f[1] = ddx12 * (sqlhf + 2*r2[1]*(vpot[1] - E))
        for i in 1:Nrmesh-1
            f[i+1] = ddx12 * (sqlhf + 2*r2[i+1]*(vpot[i+1] - E)) # Ha
            if f[i+1] == 0.0
                f[i+1] = 1.0e-20
            end
            if f[i+1] != copysign(f[i+1], f[i])
                icl = i
            end
        end
        println("icl=$icl r=$(r[icl])")
        if icl < 0 || icl >= Nrmesh - 2
            eup = E
            E = 0.5 * (eup + elw)
            continue
        end

        f .= 1.0 .- f
        fill!(y, 0.0)

        # Outward integration from origin
        y[1] = r[1]^(l+1) * (1.0 - 2.0 * Zatom * r[1] / x2l2) / sqr[1]
        y[2] = r[2]^(l+1) * (1.0 - 2.0 * Zatom * r[2] / x2l2) / sqr[2]

        ncross = 0
        for i in 1:icl-1
            y[i+2] = ((12.0 - 10.0 * f[i+1]) * y[i+1] - f[i] * y[i]) / f[i+2]
            if y[i+1] != copysign(y[i+1], y[i+2])
                ncross += 1
            end
        end
        fac = y[icl+1]

        if ncross != nodes
            if ncross > nodes
                eup = E
            else
                elw = E
            end
            E = 0.5 * (eup + elw)
            continue
        end

        # Inward integration from large r
        y[end] = dx
        y[end-1] = (12.0 - 10.0 * f[end]) * y[end] / f[end-1]

        for i in (Nrmesh-2):-1:(icl+1)
            y[i] = ((12.0 - 10.0 * f[i+1]) * y[i+1] - f[i+2] * y[i+2]) / f[i]
            # Don't let it become too big
            if y[i] > 1.0e10
                for j in Nrmesh-1:-1:(i-1)
                    y[j+1] /= y[i]
                end
            end
        end

        # Match at classical turning point
        fac = fac / y[icl+1]
        y[icl+1:end] .*= fac

        # Normalization
        norm = sum(y .* y .* r2 .* dx)
        norm = sqrt(norm)
        y ./= norm

        norm = sum(y .* y .* r2 .* dx)
        println("After normalization: norm = ", norm)

        # Cusp correction at matching point
        i = icl
        ycusp = (y[i] * f[i] + f[i+2] * y[i+2] + 10.0 * f[i+1] * y[i+1]) / 12.0
        dfcusp = f[i+1] * (y[i+1] / ycusp - 1.0)

        de = 0.5*dfcusp / ddx12 * ycusp * ycusp * dx # Ha
        if de > 0.0
            elw = E
        end
        if de < 0.0
            eup = E
        end
        E_old = E

        E = max(min(E + de, eup), elw)

        @printf("de = %18.10e\n", de)
        println("diff E = ", abs(E-E_old))
        if abs(de) < TOL
            converged = true
            break
        end
    end

    if !converged
        if ncross != nodes
            println(E, " ", elw, " ", eup, " ", ncross, " ", nodes, " ", icl)
        else
            println(E, " ", de)
        end
        error("error in solve_sch_rad: too many iterations")
    else
        @printf(" convergence achieved at iter #%3d de = %10.4e\n", iter, de)
    end

    return E, y
end


function solve_sch_rad_L(
    n, l, Nrmesh, dx, r, sqr, r2, vpot, Zatom;
    maxiter = 200,
    TOL = 1.0e-10
)

    f = zeros(Float64, Nrmesh)
    y_out = zeros(Float64, Nrmesh)
    y_in = zeros(Float64, Nrmesh)

    ddx12 = dx*dx/12.0
    sqlhf = (l + 0.5)^2
    x2l2  = 2*l + 2

    # Energy bounds
    eup = vpot[end]
    elw = minimum( sqlhf ./ (2.0 .* r2) .+ vpot )
    if eup - elw < TOL
        error("solve_sch_rad_L: lower and upper bounds are equal")
    end

    E = 0.5 * (elw + eup)
    converged = false
    iter = 0

    nodes = n - l - 1
    ncross = -1
    icl = -1

    for kkk in 1:maxiter
        iter = kkk

        # Set up f-function and find classical turning point
        icl = -1
        f[1] = ddx12 * (sqlhf + 2*r2[1]*(vpot[1] - E))
        for i in 1:Nrmesh-1
            f[i+1] = ddx12 * (sqlhf + 2*r2[i+1]*(vpot[i+1] - E))
            if f[i+1] == 0.0
                f[i+1] = 1.0e-20
            end
            if f[i+1] != copysign(f[i+1], f[i])
                icl = i
            end
        end
        println("icl=$icl r=$(r[icl])")

        if icl < 0 || icl >= Nrmesh - 2
            eup = E
            E = 0.5 * (eup + elw)
            continue
        end

        f .= 1.0 .- f
        fill!(y_out, 0.0)
        fill!(y_in, 0.0)

        # Outward integration from origin
        y_out[1] = r[1]^(l+1) * (1.0 - 2.0 * Zatom * r[1] / x2l2) / sqr[1]
        y_out[2] = r[2]^(l+1) * (1.0 - 2.0 * Zatom * r[2] / x2l2) / sqr[2]
        ncross = 0
        for i in 1:icl-1
            y_out[i+2] = ((12.0 - 10.0 * f[i+1]) * y_out[i+1] - f[i] * y_out[i]) / f[i+2]
            if y_out[i+1] != copysign(y_out[i+1], y_out[i+2])
                ncross += 1
            end
        end

        # Node-count check
        if ncross != nodes
            if ncross > nodes
                eup = E
            else
                elw = E
            end
            E = 0.5 * (eup + elw)
            continue
        end

        println("Found needed nodes: eup=$eup elw=$elw E=$E")

        # Inward integration from large r
        y_in[end] = dx
        y_in[end-1] = (12.0 - 10.0 * f[end]) * y_in[end] / f[end-1]

        for i in (Nrmesh-2):-1:(icl+1)
            y_in[i] = ((12.0 - 10.0 * f[i+1]) * y_in[i+1] - f[i+2] * y_in[i+2]) / f[i]

            # Rescale **including** y_in[i] itself
            if abs(y_in[i]) > 1.0e10
                scale = 1.0 / y_in[i]
                for j in i:Nrmesh
                    y_in[j] *= scale
                end
            end
        end

        # Matching point
        m = icl + 1
        if m == 1 || m == Nrmesh
            error("Matching point at boundary")
        end

        # This is not really necessary
        #ss = y_out[m]/y_in[m]
        ## scale only m+1 and m
        #y_in[m] *= ss
        #y_in[m+1] *= ss
        #println("y_out[m] = ", y_out[m], " y_in[m] = ", y_in[m])
        
        # Logarithmic derivative mismatch L(E)
        d_out = (y_out[m] - y_out[m-1]) / (r[m] - r[m-1])
        d_in  = (y_in[m+1] - y_in[m]) / (r[m+1] - r[m])

        if abs(y_out[m]) < 1e-30 || abs(y_in[m]) < 1e-30
            error("Zero value at matching point")
        end

        L = d_out/y_out[m] - d_in/y_in[m]

        # Correct sign: L > 0 ⇒ energy is too low ⇒ raise lower bound
        if L > 0.0
            elw = E
        else
            eup = E
        end

        e_new = 0.5 * (elw + eup)

        diff_E = abs(eup - elw)
        @printf("iter %3d  E = %18.10e  L(E) = %12.4e diff_E = %12.4e\n", kkk, e_new, L, diff_E)
        #println("Derivative mismatch: d_out=$d_out d_in=$d_in diff=$(abs(d_out-d_in))")

        if abs(L) < TOL || diff_E < TOL
            E = e_new
            converged = true
            break
        end

        E = e_new
    end

    if !converged
        error("solve_sch_rad_L: too many iterations")
    end

    # Construct matched wavefunction (optional)
    m = icl + 1
    fac = y_out[m] / y_in[m]
    y_in[m:end] .*= fac
    y = copy(y_out)
    y[m:end] = y_in[m:end]

    # Normalization
    norm = sum(y .* y .* r2 .* dx)
    norm = sqrt(norm)
    y ./= norm
    norm = sum(y .* y .* r2 .* dx)
    println("After normalization: norm = ", norm)

    # Calculate mismatch again
    m = icl + 1
    if m == 1 || m == Nrmesh
        error("Matching point at boundary")
    end    
    # Logarithmic derivative mismatch L(E)
    d_out = (y[m] - y[m-1]) / (r[m] - r[m-1])
    d_in  = (y[m+1] - y[m]) / (r[m+1] - r[m])
    L = d_out/y[m] - d_in/y[m]
    println("At the end: L = $L d_in=$d_in d_out=$d_out")
    #
    return E, y
end

function main_debug()
    Zatom = 4.0
    n = 4
    l = 1
    @assert Zatom >= 1.0
    @assert n >= 1
    @assert n <= Int(Zatom)
    @assert l < n
    println("Zatom = $Zatom n = $n l = $l")

    zmesh = Zatom
    rmax = 100.0
    xmin = -8.0
    dx = 0.01

    Nrmesh = Int(floor((log(zmesh * rmax) - xmin) / dx))
    r, sqr, r2 = initialize_mesh(zmesh, xmin, dx, Nrmesh)
    vpot = initialize_V_atom(Zatom, r, Nrmesh)

    E_nl1, y1 = solve_sch_rad(n, l, Nrmesh, dx, r, sqr, r2, vpot, Zatom)
    E_nl2, y2 = solve_sch_rad_L(n, l, Nrmesh, dx, r, sqr, r2, vpot, Zatom)

    @printf("Energies in Ha\n")
    @printf("E_nl1=%15.8f E_nl2=%15.8f,  Analytic=%15.8f\n", E_nl1, E_nl2, -0.5*(Zatom/n)^2)

    fig =  Plots.plot(r, y1)
    Plots.plot!(r, y2)
    Plots.xlims!(0.0, 10.0)
    @exfiltrate
    return fig
end

#main()