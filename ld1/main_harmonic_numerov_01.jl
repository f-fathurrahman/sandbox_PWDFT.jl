using Printf

xmax = 10.0
Nxmesh = 500

x = zeros(Float64, Nxmesh)
y = zeros(Float64, Nxmesh)
p = zeros(Float64, Nxmesh)
Vpot = zeros(Float64, Nxmesh)
f = zeros(Float64, Nxmesh)

dx = xmax/(Nxmesh-1) 
ddx12 = dx*dx/12.0

# set up the potential (must be even w.r.t. x=0)
for i in 1:Nxmesh
    x[i] = (i-1)*dx
    Vpot[i] = 0.5*x[i]*x[i]
end

Nnodes = 2
@assert Nnodes >= 0

#
# set initial lower and upper bounds to the eigenvalue
#
eup = maximum(Vpot)
elw = minimum(Vpot)

# search eigenvalues with bisection (max 1000 iterations)
e = 0.5 * (elw + eup)
n_iter = 1000


icl = -1 # put here so that it can be accessed outside the loop

#for kkk in 1:n_iter
for kkk in 1:50
    #
    # set up the f-function used by the Numerov algorithm
    # and determine the position of its last crossing, i.e. change of sign
    # f < 0 means classically allowed   region
    # f > 0 means classically forbidden region
    #
    # classically forbidden region (for which V(x) > E)
    #
    f[1] = ddx12*( 2.0*( Vpot[1] - e) )
    #
    # TODO: make this a function: find_icl(...)
    # We only need array Vpot - E
    for i in 2:Nxmesh
        f[i] = ddx12*2.0*( Vpot[i] - e )
        # beware: if f(i) is exactly zero the change of sign is not observed
        # the following line is a trick to prevent missing a change of sign 
        # in this unlikely but not impossible case:
        # XXX: what's this?
        if f[i] == 0.0
            f[i] = 1.e-20
        end
        # store the index 'icl' where the last change of sign has been found
        if sign(f[i]*f[i-1]) < 0
            icl = i
        end
    end

    if icl >= Nxmesh-2
        println("Error: last change of sign too far")
    elseif icl < 1
        println("Error: no classical turning point")
    end

    #println("icl = ", icl)


    #
    # f(x) as required by the Numerov algorithm
    #
    for i in 1:Nxmesh
        f[i] = 1 + 2*( e - Vpot[i] )*ddx12
    end
    # The wave function, reset to zeros
    fill!(y, 0.0)
    
    # determination of the wave-function in the first two points 
    hnodes = floor(Int64, Nnodes/2) # half nodes
    # 1/2 -> 0

    #
    # beware the integer division: 1/2 = 0 !
    # if nodes is even, there are 2*hnodes nodes
    # if nodes is odd,  there are 2*hnodes+1 nodes (one is in x=0)
    # hnodes is thus the number of nodes in the x > 0 semi-axis (x=0 is excepted)


    if 2*hnodes == Nnodes
        # even number of nodes: wavefunction is even
        y[1] = 1.0
        # assume f[-1] = f[1] (XXX FIXME)
        y[2] = ( 12.0 - 10.0*f[1] )*y[1]/(2*f[2])
    else
        # odd number of nodes: wavefunction is odd
        y[1] = 0.0
        y[2] = dx
    end

    # We are doing outward integration (because we only consider half domain)
    # outward integration and count number of crossings
    ncross = 0
    for i in 2:Nxmesh-1
        y[i+1] = ( (12.0 - 10.0*f[i])*y[i] - f[i-1]*y[i-1] )/ f[i+1]
        # Check for changed sign
        if sign(y[i] * y[i+1]) < 0
            ncross = ncross + 1
        end
    end

    @printf("%3d %18.10f %3d %3d\n", kkk, e, ncross, hnodes)

    # if iterating on energy: check number of crossings
    if ncross > hnodes
        # Too many crossings: current energy is too high
        # lower the upper bound
        eup = e
        # eup is set to e so that e cannot change to larger than eup
    else 
        # Too few (or correct) number of crossings:
        # current energy is too low, raise the lower bound
        elw = e
    end
    # XXX we might want to stop based on number of crossings (?)

    # New trial value: (using bisection)
    e = 0.5 * (eup + elw)
    # Convergence criterion:
    Δe = abs(eup - elw)
    if Δe < 1.e-10
        println("Converged because Δe = ", Δe)
        break
    end

end


#
# ---- convergence has been achieved (or it wasn't required) -----
# Note that the wavefunction is not normalized: 
# the problem is the divergence at large |x| 
#
# Calculation of the classical probability density for energy e:
#
nrm = 0.0
@views p[icl:end] .= 0.0
for i in 1:icl
    arg = (e - x[i]^2/2.0)
    if arg > 0.0
        p[i] = 1.0/sqrt(arg)
    else
        p[i] = 0.0
    end
    nrm += 2.0 * dx * p[i]
end
# The point at x=0 must be counted once:
nrm -= dx*p[1] # fix double counting

# Normalize p(x) so that  Int p(x)dx = 1
@views p[1:icl-1] .*= (1/nrm)


