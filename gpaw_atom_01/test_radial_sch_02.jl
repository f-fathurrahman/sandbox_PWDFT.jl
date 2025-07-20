using Infiltrator

using LinearAlgebra 
using Printf

using Plots, PlotThemes
theme(:dark)

using PWDFT: LibxcXCCalculator, calc_epsxc_Vxc_LDA!

includet("AERadialGrid.jl")
includet("my_functions.jl")

function main_radial_02()

    atsymb = "Si"
    Z = 14

    # Input configuration
    n_j =  [1, 2, 2, 3, 3] # main/principal quantum number
    l_j =  [0, 0, 1, 0, 1] # orbital
    f_j =  [2, 2, 6, 2, 2] # occupation
    e_j =  [-65.184426, -5.075056, -3.514938, -0.398139, -0.153293]

    maxnodes = maximum([n - l - 1 for (n, l) in zip(n_j, l_j)])
    #print("maxnodes = ", maxnodes)

    gpernode = 150 # default?
    β = 0.4 # parameter for radial grid
    NradialPoints = gpernode*(maxnodes + 1)

    a = β/NradialPoints
    b = 1/NradialPoints
    ae_grid = AERadialGrid(a, b, npts=NradialPoints)

    r = ae_grid.r
    dr = ae_grid.dr
    d2gdr2 = ae_grid.d2gdr2

    Nstates = length(n_j)
    # Radial wave functions multiplied by radius:
    u_j = Vector{Vector{Float64}}(undef, Nstates)
    for ist in 1:Nstates
        u_j[ist] = zeros(Float64, NradialPoints)
    end

    # Effective potential multiplied by radius:
    vr = zeros(Float64, NradialPoints)

    # Electron density:
    rhoe = zeros(Float64, NradialPoints)

    vHr = zeros(Float64, NradialPoints)
    vXC = zeros(Float64, NradialPoints)
    epsxc = zeros(Float64, NradialPoints)

    # Initialize starting wavefunctions and calculate density from them
    init_wavefunc!(atsymb, r, dr, l_j, e_j, u_j)
    calc_density!(r, f_j, u_j, rhoe)

    solve_radial_hartree!(0, rhoe .* r .* dr, r, vHr)

    # add potential from nuclear point charge (v = -Z / r)
    vHr .-= Z # vHr is vH times r, so the potential to added is (-Z/r)*r => -Z

    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)

    fill!(vXC, 0.0)
    calc_epsxc_Vxc_LDA!(xc_calc, rhoe, epsxc, vXC)

    Exc = radial_integrate(ae_grid, epsxc .* rhoe)
    @. vr = (vHr + vXC*r)

    for ist in 1:Nstates
        fill!(u_j[ist], 0.0)
    end

    scalarrel = true

    c2 = @. -(r/dr)^2 # beware this can become a matrix without @.
    c10 = @. -d2gdr2 * r^2  # first part of c1 vector

    if scalarrel
        r2dvdr = zeros(Float64, NradialPoints)
        radial_derivative!(ae_grid, vr, r2dvdr) # calc dvdr
        r2dvdr[:] .*= r
        r2dvdr[:] .-= vr
    else
        r2dvdr = nothing
    end
    #println("sum abs r2dvdr = ", sum(abs.(r2dvdr)))
    # XXX: r2dvdr is local to this function

    # solve for each quantum state separately
    Nstates = length(n_j)
    ist = 4
    
    println("\nSolving radial Schroedinger equation for ist = $ist")

    n = n_j[ist]
    l = l_j[ist]
    E = -0.5573946000000001
    u = u_j[ist]
    #
    # initial radial integration
    nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)


    @exfiltrate
end
