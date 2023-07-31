using Printf
using Plots, PlotThemes

theme(:dark)

function find_icl(E::Float64, V::Vector{Float64})

    N = size(V, 1)

    fm1 = V[1] - E
    if fm1 == 0.0
        fm1 = 1e-20
    end

    icl = -1
    for i in 2:N
        f = V[i] - E
        #XXX Do let let f become too small (?)
        if f == 0.0
            f = 1e-20
        end
        if sign(f*fm1) < 0
            icl = i
        end
        fm1 = f
    end

    return icl
end


function integ_numerov_sch_inward!(
    icl::Int64,
    E::Float64,
    V::Vector{Float64},
    Δx::Float64,
    y::Vector{Float64};
    parity=:odd
)



    Nmesh = size(V,1)
    @assert Nmesh == size(y,1)

    # XXX: Need to calculate beyond icl ?
    f = zeros(Float64, Nmesh)
    for i in 1:Nmesh
        f[i] = 1 + 2*( E - V[i] )*Δx^2/12
    end

    y[Nmesh] = Δx
    y[Nmesh-1] = ( 12.0 - 10.0*f[Nmesh] )*y[Nmesh]/f[Nmesh-1]
    nrm_guard = 1e100
    for i in range(Nmesh-1,icl+1,-1)
        y[i-1] = ( (12.0 - 10.0*f[i] )*y[i] - f[i+1]*y[i+1])/f[i-1]
        # the following lines prevent overflows if starting from too far
        if abs(y[i-1]) > nrm_guard
            @views y[i-1:Nmesh] .= y[i-1:Nmesh] / nrm_guard
        end
    end

    return

end






function integ_numerov_sch_outward!(
    icl::Int64,
    E::Float64,
    V::Vector{Float64},
    Δx::Float64,
    y::Vector{Float64};
    parity=:odd
)


    Nmesh = size(V,1)
    @assert Nmesh == size(y,1)

    # XXX: Need to calculate beyond icl ?
    f = zeros(Float64, Nmesh)
    for i in 1:Nmesh
        f[i] = 1 + 2*( E - V[i] )*Δx^2/12
    end

    # Two possible starting points
    # Is this important (?)
    if parity == :odd
        y[1] = 1.0
        y[2] = ( 12.0 - 10.0*f[1] )*y[1]/(2*f[2])
    elseif parity == :even
        y[1] = 0.0
        y[2] = Δx
    else
        error("Unknown parity: ", parity)
    end

    Ncross = 0
    for i in 2:icl-1
        y[i+1] = ( (12.0 - 10.0*f[i])*y[i] - f[i-1]*y[i-1] )/ f[i+1]
        # Check for changed sign
        if sign(y[i] * y[i+1]) < 0
            println("Found crossing: ", i, "--", i + 1)
            Ncross = Ncross + 1
        end
    end

    return Ncross
end


function main()

    a = -5.0
    b =  5.0
    Npoints = 500
    Δx = (b - a)/(Npoints-1)
    xgrid = zeros(Float64, Npoints)
    for i in 1:Npoints
        xgrid[i] = a + (i-1)*Δx
    end
    println("xgrid[1] = ", xgrid[1])
    println("xgrid[end] = ", xgrid[end])

    Vpot = zeros(Float64, Npoints)
    for i in 1:Npoints
        Vpot[i] = 0.5*xgrid[i]^2
    end

    #plot(xgrid, Vpot, label="Potential", linewidth=2, dpi=150)
    #title!("Harmonic Potential")
    #savefig("IMG_pot.png")

    E = 0.1
    psi = zeros(Float64, Npoints)
    icl = find_icl(E, Vpot)
    println("icl = ", icl)
    #Ncross = integ_numerov_sch!( E, Vpot, Δx, psi, parity=:even )
    Ncross = integ_numerov_sch_outward!( icl, E, Vpot, Δx, psi, parity=:odd )
    println("Ncross = ", Ncross)

    #nrm = sum(psi .* psi) * Δx
    #psi *= (1/sqrt(nrm))
    #nrm = sum(psi .* psi) * Δx
    #println("nrm psi = ", nrm)

    plot(xgrid[1:icl], psi[1:icl], label="psi1", linewidth=2, dpi=200)
    savefig("IMG_psi1.png")

end

main()

