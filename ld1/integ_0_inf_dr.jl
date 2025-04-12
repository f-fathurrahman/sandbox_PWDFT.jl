function integ_0_inf_dr(f, grid, Nrmesh, nst)
    # XXX Need Nrmesh?
    #=
    integral of f from 0 to infinity
    f is given on a logarithmic mesh. 
    f(r) is assumed to be proportional to r**nst for small r
    =#

    @assert Nrmesh <= grid.Nrmesh

    fs = zeros(Float64, 4)
    b = zeros(Float64, 4)
    # series development: contribution for small r
    for i in 1:4
        fs[i] = f[i] / grid.r[i]^nst
    end
    radial_grid_series!( fs, grid.r, grid.r2, b )
    
    res = ( b[1]/(nst+1) + grid.r[1]*( b[2]/(nst+2) + grid.r[1]*b[3]/(nst+3)) ) * grid.r[1]^(nst+1)

    # simpson integration (logarithmic mesh: dr ==> r dx)
    #
    ss = 0.0
    for i in range(1, stop=Nmesh-2, step=2)
        ss1 += f[i]*grid.r[i] + 4.0*f[i+1]*grid.r[i+1] + f[i+2]*grid.r[i+2]
    end
    res += ss1*grid.dx/3.0
    return res
end

