# Initial minimization
function minimize_electrons!( Ham, psiks; NiterMax=200, etot_conv_thr=1e-8 )
    KS_solve_Emin_PCG!( Ham, psiks,
        skip_initial_diag=true, etot_conv_thr=etot_conv_thr,
        NiterMax=NiterMax
    )
    forces = calc_forces( Ham, psiks )
    return sum(Ham.energies), forces
end

# μ is fictitious electron mass
function calc_elec_energies!(Ham, psiks, μ, dpsiks)
    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, psiks)
    energies = calc_energies(Ham, psiks)
    Etot = sum(energies)
    Ekin_elec = 0.5*μ*dot(dpsiks,dpsiks)
    return Etot, Ekin_elec
end 

# RATTLE algorithm
function calc_X_matrix(Ctilde::Array{ComplexF64,2}, C::Array{ComplexF64,2})
    A = Ctilde' * Ctilde
    B = C' * Ctilde
    X = 0.5*(I - A)
    Xnew = similar(X)
    for iter in 1:100
        Xnew[:,:] = 0.5*( I - A + X*(I - B) + (I - B)*X - X*X )
        ΔX = norm(X-Xnew)
        println("norm X-Xnew = ", ΔX)
        if ΔX < 1e-10
            println("Find X: converged in iter: ", iter)
            return Xnew
        end
        @views X[:,:] = Xnew[:,:]
    end
    error("ERROR: X is not converged")
    return X
end

# Need this?
struct CarParrinelloMDParameters
    dt_fs::Float64
    dt::Float64
    μ::Float64
end