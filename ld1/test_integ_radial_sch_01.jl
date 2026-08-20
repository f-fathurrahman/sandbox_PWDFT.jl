using Printf
import LinearAlgebra

# Use import to avoid potential name clash
import Plots, PlotThemes
Plots.theme(:dark)

function find_idx_r_cl(E::Float64, V::Vector{Float64})
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
    return max(icl, 50) # at least icl should be 50, don't let it too small
end


# Effective potential for given l
function Veff(r, l, Z)
    return -Z/r + l*(l+1)/(2*r^2)
end

# Right-hand side of the first-order system: dy/dr = F(r, y)
function F(r, y, E, l, Z)
    y1, y2 = y
    factor = 2*(Veff(r, l, Z) - E)
    return [y2, factor * y1]
end

# ---- Integration routines ----
function euler_step(r, y, h, E, l, Z)
    f = F(r, y, E, l, Z)
    y_new = y + h * f
    return y_new
end

function rk4_step(r, y, h, E, l, Z)
    k1 = h * F(r,        y,        E, l, Z)
    k2 = h * F(r + h/2,  y + k1/2, E, l, Z)
    k3 = h * F(r + h/2,  y + k2/2, E, l, Z)
    k4 = h * F(r + h,    y + k3,   E, l, Z)
    y_new = y + (k1 + 2k2 + 2k3 + k4) / 6
    return y_new
end

# Integrate from r0 to rmax using chosen method
function outward_sch_integ!(rgrid, E, l, work_psi, idx_r_cl, Z)
    r0 = rgrid[1]
    # initial conditions
    work_psi[1] = r0^(l+1)
    work_psi[2] = (l+1) * r0^l
    #
    dr = rgrid[2] - rgrid[1]
    for ir in 1:idx_r_cl-1
        work_psi[1,ir+1], work_psi[2,ir+1] = rk4_step(rgrid[ir], work_psi[:,ir], dr, E, l, Z)
    end
    #println("Last point updated = ", idx_r_cl-1+1)
    return
end

function inward_sch_integ!(rgrid, E, l, work_psi, idx_r_cl, Z)
    Nrmesh = size(rgrid, 1)
    # initial conditions
    rmax = rgrid[Nrmesh]
    work_psi[1,Nrmesh] = 1.0
    krmax = sqrt(Veff(rmax, l, Z) - E)
    #println("krmax = ", krmax)
    work_psi[2,Nrmesh] = -krmax*work_psi[1,Nrmesh]
    #
    dr = rgrid[2] - rgrid[1]
    psi1_match = work_psi[1,idx_r_cl]
    psi2_match = work_psi[2,idx_r_cl]
    for ir in range(Nrmesh, stop=idx_r_cl+1, step = -1)
        work_psi[1,ir-1], work_psi[2,ir-1] = rk4_step(rgrid[ir], work_psi[:,ir], -dr, E, l, Z)
    end
    #println("Last point updated = ", idx_r_cl+1-1)
    # simply scale
    scale1 = psi1_match/work_psi[1,idx_r_cl]
    scale2 = psi2_match/work_psi[2,idx_r_cl]
    work_psi[1,idx_r_cl:end] .*= scale1
    work_psi[2,idx_r_cl:end] .*= scale2
    return
end

function debug_01()
    # initial conditions:
    # u ∝ r^{l+1}
    # u' = (l+1) r^l
    r0 = 1e-10
    rmax = 50.0
    h = 0.01
    Z = 2
    ℓ = 1
    n = 1
    rgrid = collect( range(r0, stop=rmax, step=h) )
    #E = -0.5*Z^2/n^2
    E = -0.15
    println("E = ", E)
    Nrmesh = size(rgrid, 1)

#=
    idx_r_cl = 0
    for ir in 1:Nrmesh
        if Veff(rgrid[ir], ℓ, Z) > E
            idx_r_cl = ir
            break
        end
    end
    idx_r_cl = max(50, idx_r_cl)
=#

    Vgrid = Veff.(rgrid, ℓ, Z)
    idx_r_cl = find_idx_r_cl(E, Vgrid)

    println("idx_r_cl = ", idx_r_cl)

    # work arrays
    work_psi = zeros(2, Nrmesh)
    outward_sch_integ!(rgrid, E, ℓ, work_psi, idx_r_cl, Z)
    inward_sch_integ!(rgrid, E, ℓ, work_psi, idx_r_cl, Z)

    Plots.plot(rgrid, work_psi[1,:])
    Plots.xlims!(0.0, 15.0)

    #@infiltrate
    #@exfiltrate

end