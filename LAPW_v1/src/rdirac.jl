function rdirac!(
    n::Int64, l::Int64, k::Int64,
    r, vr, evals, g0, f0;
    sol=137.035999084, max_iter=2000, tol=1e-12
)

    # ! arguments
    # real(8), intent(in) :: sol
    # integer, intent(in) :: n,l,k,nr
    # real(8), intent(in) :: r(nr),vr(nr)
    # real(8), intent(inout) :: evals
    # real(8), intent(out) :: g0(nr),f0(nr)

    nr = size(r,1)
    @assert size(vr,1) >= nr  # vr can be larger than r

    @assert k > 0
    @assert nr >= 4

    if k > n
        error("Error incompatible n=$n and k=$k")
    end

    if (k == n) && (l != k-1)
        error("Incompatible n=$n, k=$k and l=$l")
    end

    if k == l
        kpa = k
    elseif k == (l+1)
        kpa = -k
    else
        println("Error: incompatible l=$l and k=$l")
    end

    # automatic arrays
    g1 = zeros(Float64, nr)
    f1 = zeros(Float64, nr)
    fr = zeros(Float64, nr)

    de = 1.0
    nndp = 0
    iiter = 0
    #
    while true
        #
        iiter = iiter + 1
        if iiter > max_iter
            break
        end 
        # integrate the Dirac equation
        #println("calling rdiracint")
        #println("kpa = ", kpa)
        #println("nr  = ", nr)
        nn = rdiracint!(kpa, evals, r, vr, g0, g1, f0, f1, sol=sol)
        # check the number of nodes
        nnd = nn - (n-l-1)
        if nnd > 0
            evals = evals - de
        else
            evals = evals + de
        end
        #
        if iiter > 1
            if ( (nnd != 0) || ( nndp != 0) )
                if nnd*nndp <= 0
                    de = de*0.5
                else
                    de = de*1.1
                end
            end
        end
        #
        nndp = nnd
        if de < tol*(abs(evals) + 1.0)
            break
        end
    end

    if iiter > max_iter
        println("Warning(rdirac): maximum iterations exceeded")
    end

    # find effective infinity and set wavefunction to zero after that point
    # major component
    irm = nr
    for ir in 2:nr
        if (g0[ir-1]*g0[ir] < 0.0) || (g1[ir-1]*g1[ir] < 0.0)
            irm = ir
        end
    end
    g0[irm:nr] .= 0.0
    # minor component
    irm = nr
    for ir in 2:nr
        if (f0[ir-1]*f0[ir] < 0.0) || (f1[ir-1]*f1[ir] < 0.0)
            irm = ir
        end
    end
    f0[irm:nr] .= 0.0
    # normalise
    for ir in 1:nr
      fr[ir] = g0[ir]^2 + f0[ir]^2
    end

    t1 = splint(nr, r, fr)
    t1 = sqrt(abs(t1))

    if t1 > 0.0
      t1 = 1.0/t1
    else
      error("Error(rdirac): zero wavefunction")
    end
    for ir in 1:nr
        g0[ir] = t1*g0[ir]
        f0[ir] = t1*f0[ir]
    end
    return evals
end
