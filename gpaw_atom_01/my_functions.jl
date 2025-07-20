function radial_integrate(rgrid, a_xg; n=0)
    @assert n >= -2
    r = rgrid.r
    dr = rgrid.dr
    # exclude the first radial point?
    return dot( a_xg[2:end], (r.^(2 + n) .* dr)[2:end] ) * (4π)
end

# Initialize with Slater function
function init_wavefunc!(atsymb, r, dr, l_j, E_j, u_j)
    #
    Nstates = length(l_j)
    special_case = ["Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au"]
    NradialPoints = length(r)
    #
    for ist in 1:Nstates
        l = l_j[ist]
        E = E_j[ist]
        u = u_j[ist]
        if atsymb in special_case
            a = sqrt(-4.0*E)
        else
            a = sqrt(-2.0*E)
        end
        for i in 1:NradialPoints
            u[i] = r[i]^(1 + l) * exp(-a * r[i])
        end
        norm_u = dot(u.^2, dr) # use radial_wavefunc_norm
        u[:] .*= 1.0 / sqrt(norm_u)
    end
    return
end

# Return the electron charge density divided by 4 pi
function calc_density!(r, f_j, u_j, rhoe)
    NradialPoints = length(rhoe)
    Nstates = length(f_j)
    pref = 1/(4π)
    fill!(rhoe, 0.0)
    for ist in 1:Nstates, i in 2:NradialPoints
        if abs(u_j[ist][i]) >= 1e-160
            rhoe[i] += pref * f_j[ist] * u_j[ist][i]^2 / r[i]^2
        end
    end
    rhoe[1] = rhoe[2]
    return
end

# input n*r*dr, output vr=vh*r
function solve_radial_hartree!(
    l::Int64, nrdr, r, vr
)

    NradialPoints = length(r)    
    p = 0.0
    q = 0.0
    fill!(vr, 0.0)
    for i in range(NradialPoints, stop=1, step=-1)
        rl = r[i]^l
        dp = nrdr[i]/rl
        rlp1 = rl * r[i]
        dq = nrdr[i] * rlp1
        vr[i] = (p + 0.5*dp)*rlp1 - (q + 0.5*dq)/rl
        p += dp
        q += dq
    end

    #println("p = ", p, " q = ", q)
    
    vr[1] = 0.0
    f = 4π / (2*l + 1)
    #for (int g = 1; g < M; g++)
    for i in 2:NradialPoints
        vr[i] = f * (vr[i] + q / r[i]^l)
    end

    return
end

# Finite-difference derivative of radial function.
function radial_derivative!(rgrid, n_g, dndr_g)
    #
    npts = rgrid.npts
    #
    dndr_g[1] = n_g[2] - n_g[1] # forward difference
    #
    # centered difference
    for i in 2:npts-1
        dndr_g[i] = 0.5*(n_g[i+1] - n_g[i-1])
    end
    dndr_g[npts] = n_g[npts] - n_g[npts-1] # backward difference
    @. dndr_g /= rgrid.dr
    return
end

function radial_wavefunc_norm(u, dr)
    idx = abs.(u) .>= 1e-160
    return dot( u[idx].^2, dr[idx] ) 
end

includet("solve_radial_sch.jl")
includet("integ_radial_sch.jl")