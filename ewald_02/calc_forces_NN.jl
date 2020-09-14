using SpecialFunctions: erfc
import PWDFT: calc_forces_NN, calc_forces_NN!

function calc_forces_NN( pw::PWGrid, atoms::Atoms )
    F_NN = zeros(3, atoms.Natoms)
    calc_forces_NN!( atoms, F_NN )
    return F_NN
end

function calc_forces_NN!( pw::PWGrid, atoms::Atoms, F_NN::Array{Float64,2} )
    return calc_forces_NN!( pw, atoms, atoms.Zvals, F_NN )
end

function _Herfc(x)
    return -2*exp(-x^2)/sqrt(pi) - erfc(x)/x
end

function calc_forces_NN!(
    pw::PWGrid,
    atoms::Atoms,
    Zvals::Array{Float64,1},
    F_NN::Array{Float64,2}
)

    LatVecs = pw.LatVecs  
    Ω = pw.CellVolume
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # Atomic positions
    tau = atoms.positions

    charge = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        charge = charge + Zvals[isp]
    end

    G2_max = pw.ecutrho
    α = 1.0
    upperbound = charge^2 * sqrt(α/π) * erfc(sqrt(G2_max/4.0/α))
    while upperbound > 1e-6
        α = α - 0.1
        if alpha <= 0.0
            error("Optimal α is not found")
        end
        upperbound = charge^2 * sqrt(α/π) * erfc(sqrt(G2_max/4.0/α))
    end

    println("α = ", α)

    Ng = pw.gvec.Ng
    G = pw.gvec.G
    G2 = pw.gvec.G2
    aux = zeros(ComplexF64,Ng)
    for ia in 1:Natoms
        isp = atm2species[ia]
        for ig in 2:Ng
            GX = tau[1,ia]*G[1,ig] + tau[2,ia]*G[2,ig] + tau[3,ia]*G[3,ig]
            Sf = cos(GX) + im*sin(GX) # conj
            aux[ig] = aux[ig] + Zvals[isp] * Sf
        end
    end

    println("Ng = ", pw.gvec.Ng)
    println("sum aux = ", sum(aux))
    println("some aux")
    for ig in 2:10
        @printf("%8d %18.10f %18.10f %18.10f\n", ig, G2[ig], real(aux[ig]), imag(aux[ig]))
    end

    # Treat 2d cutoff is skipped

    # skip G2=0
    for ig in 2:Ng
        aux[ig] = aux[ig] * exp(-G2[ig]/α/4.0) / G2[ig]
    end

    println("some aux after modified")
    for ig in 2:10
        @printf("%8d %18.10f %18.10f %18.10f\n", ig, G2[ig], real(aux[ig]), imag(aux[ig]))
    end

    println("Ω = ", Ω)
    alat = 16.0 # HARDCODED
    Gfact = 2.0

    F_NN_G = zeros(3,Natoms)

    for ia in 1:Natoms
        isp = atm2species[ia]
        for ig in 2:Ng # from 1?
            GX = G[1,ig]*tau[1,ia] + G[2,ig]*tau[2,ia] + G[3,ig]*tau[3,ia]
            sumnb = cos(GX)*imag(aux[ig]) - sin(GX)*real(aux[ig])
            F_NN_G[1,ia] = F_NN_G[1,ia] + G[1,ig] * sumnb #/ (2π/alat)
            F_NN_G[2,ia] = F_NN_G[2,ia] + G[2,ig] * sumnb #/ (2π/alat)
            F_NN_G[3,ia] = F_NN_G[3,ia] + G[3,ig] * sumnb #/ (2π/alat)
        end
        #fact = -Zvals[isp] * Gfact / pw.CellVolume
        #fact = -Zvals[isp] * Gfact * (2*pi)^2 / pw.CellVolume / alat
        fact = -4π * Zvals[isp]  / pw.CellVolume
        for i in 1:3
            F_NN_G[i,ia] = F_NN_G[i,ia]*fact
        end
    end

    fact = -Zvals[1] * Gfact / pw.CellVolume
    @printf("fact = %18.10e\n", fact)

    alat = 16.0
    fact = -Zvals[1] * Gfact * 2.0 * (2*pi)^2 / pw.CellVolume / alat
    @printf("fact = %18.10e\n", fact)

    println("Reciprocal space sum contribution: ")
    display(F_NN_G'); println()

    F_NN[:] = F_NN_G[:] #+ F_NN_R[:]

    return
end