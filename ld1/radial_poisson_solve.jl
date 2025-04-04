# Solution of the Poisson's equation on a radial (logarithmic) grid
function radial_poisson_solve!(k, nst, grid, f, vh)
    # k,  ! input: the k of the equation
    # nst ! input: at low r, f goes as r**nst
    #real(DP), intent(in) :: f(Nrmesh)  input: the 4\pi r2 \rho function
  
    Nrmesh = grid.Nrmesh
    d = zeros(Float64,Nrmesh)
    e = zeros(Float64,Nrmesh)
    
    # Find the series expansion of the solution close to r=0
    k21 = 2*k + 1
    nk1 = nst - k -1
    if nk1 <= 0
        error("Error in radial_solve_poison: nk1 <= 1")
    elseif nk1 >= 3
        #else if(nk1.ge.4) then
        # not sure whether the following is really correct, but the above wasn't
        c2 = 0.0
        c3 = 0.0
    else
        e[1] = 0.0
        for i in 1:4
           d[i] = -k21*f[i]/grid.r[i]^nst
        end
        idx_e = nk1:(nk1+3)
        @views radial_grid_series!(d, grid.r, grid.r2, e[idx_e])
        c2 = e[1]/(4.0*k + 6.0)
        c3 = e[2]/(6.0*k + 12.0)
    end

    # Set the main auxiliary parameters
    ch = grid.dx^2/12.0
    xkh2 = ch*(k + 0.5)^2
    ei = 1.0 - xkh2
    di = -(2.0 + 10.0*xkh2)

    # Set the diagonal and the off diagonal elements of the 
    # linear system, compute a part of the right hand side 
    for i in 2:Nrmesh
       d[i] = -di
       e[i] = -ei
       vh[i] = k21*ch*grid.sqrtr[i]*f[i]
    end

    # Use the boundary condition to eliminate the value of the 
    # solution in the first point from the first equation. This 
    # part for the diagonal element
    f1 = (grid.sqrtr[1]/grid.sqrtr[2])^k21
    d[2] = d[2] - ei*f1
    
    # Use the boundary condition to eliminate the value of the 
    # solution in the last point from the last equation
    fn = ( grid.sqrtr[Nrmesh-1]/grid.sqrtr[Nrmesh] )^k21
    d[Nrmesh-1] = d[Nrmesh-1] - ei*fn

    # In the first point vh(1) has the same definition as in the other points
    vhim1 = k21*ch*grid.sqrtr[1]*f[1]
    
    # Compute the right hand side using the auxiliary quantity vh(i).
    for i in 2:Nrmesh-1
        vhi = vh[i]
        vh[i] = vhim1 + 10.0*vhi + vh[i+1]
        vhim1 = vhi
    end
    # Use the boundary condition to eliminate the value of the solution in the 
    # first point from the first equation. This part for the right hand side.
    #
    vh[2] = vh[2] - ei*grid.sqrtr[1]^k21 * (
            c2*(grid.r2[2] - grid.r2[1]) + 
            c3*(grid.r[2]^3 - grid.r[1]^3) )

    # solve the linear system with lapack routine dptsv
    idx_start = 2
    idx_stop = 2 + (Nrmesh-2) - 1
    idx_D = idx_start:idx_stop
    idx_E = idx_start:(idx_stop-1)

    #@infiltrate

    @views LinearAlgebra.LAPACK.ptsv!( d[idx_D], e[idx_E], vh[idx_D] )
    # Original call call dptsv(mesh-2,1,d(2),e(2),vh(2),mesh-2,ierr)
    # ptsv!(D, E, B)
    #if (ierr.ne.0) call upf_error('hartree', 'error in lapack', ierr)
  
    # Set the value of the solution at the first and last point
    # First, find c0 from the solution in the second point
    c0 = vh[2]/grid.sqrtr[2]^k21 - c2*grid.r2[2] - c3*grid.r[2]*grid.r2[2]

    # and then use the series expansion at the first point
    vh[1] = grid.sqrtr[1]^k21*( c0 + c2*grid.r2[1] + c3*grid.r[1]^3 )

    # the solution at the last point is given  by the boundary 
    # condition
    vh[Nrmesh] = vh[Nrmesh-1]*fn

    # The solution must be divided by r (from the equation) 
    # and multiplied by the square root of r (from the log 
    # mesh transformation)
    for i in 1:Nrmesh
       vh[i] = vh[i] / grid.sqrtr[i]
    end

    #@infiltrate

    return
end
