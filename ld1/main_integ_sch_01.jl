using Printf
using Plots, PlotThemes

theme(:dark)

function integ_numerov_sch_outward!(
    E::Float64,
    V::Vector{Float64},
    Δx::Float64,
    y::Vector{Float64};
    LARGE=1e5,
    parity=:odd
)


    Nmesh = size(V,1)
    @assert Nmesh == size(y,1)

    f = zeros(Float64, Nmesh)
    for i in 1:Nmesh
        f[i] = 1 + 2*( E - V[i] )*Δx^2/12
    end

    # Two possible starting points
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
    iend = Nmesh
    for i in 2:Nmesh-1
        y[i+1] = ( (12.0 - 10.0*f[i])*y[i] - f[i-1]*y[i-1] )/ f[i+1]
        # Check for changed sign
        if sign(y[i] * y[i+1]) < 0
            println("Found crossing: ", i, "--", i + 1)
            Ncross = Ncross + 1
        end
        if abs(y[i+1]) > LARGE
            iend = i + 1
            break
        end
    end

    return Ncross, iend
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

    E = 1.5
    psi = zeros(Float64, Npoints)
    #Ncross = integ_numerov_sch!( E, Vpot, Δx, psi, parity=:even )
    Ncross, iend = integ_numerov_sch_outward!( E, Vpot, Δx, psi, parity=:even )
    println("Ncross = ", Ncross)
    println("iend = ", iend)

    #nrm = sum(psi .* psi) * Δx
    #psi *= (1/sqrt(nrm))
    #nrm = sum(psi .* psi) * Δx
    #println("nrm psi = ", nrm)

    plot(xgrid[1:iend], psi[1:iend], label="psi1", linewidth=2, dpi=200)
    savefig("IMG_psi1.png")

end

main()

