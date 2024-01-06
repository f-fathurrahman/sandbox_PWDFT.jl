#
# adapted from rgen.f90 of QE-6.6
#
function gen_neighbor_shells!(
    dtau, rmax::Float64, mxr::Int64,
    LatVecs, RecVecs,
    r, r2
)
    #
    # generates neighbours shells (cartesian, in units of lattice parameter)
    # with length < rmax,and returns them in order of increasing length:
    #    r(:) = i*a1(:) + j*a2(:) + k*a3(:) - dtau(:),   r2 = r^2
    # where a1, a2, a3 are primitive lattice vectors. Other input variables:
    #   mxr = maximum number of vectors
    #   at  = lattice vectors ( a1=at(:,1), a2=at(:,2), a3=at(:,3) )
    #   bg  = reciprocal lattice vectors ( b1=bg(:,1), b2=bg(:,2), b3=bg(:,3) )
    # Other output variables:
    #   nrm = the number of vectors with r^2 < rmax^2
    #
    # USE kinds, ONLY : DP
    # !
    # IMPLICIT NONE
    # INTEGER, INTENT(in) :: mxr
    # INTEGER, INTENT(out):: nrm
    # REAL(DP), INTENT(in) :: at(3,3), bg(3,3), dtau(3), rmax
    # REAL(DP), INTENT(out):: r(3,mxr), r2(mxr)

    nrm = 0
    SMALL = eps()
    if rmax <= SMALL
        return nrm
    end
  
    ds = zeros(Float64, 3)
    dtau0 = zeros(Float64, 3)

    # bring dtau into the unit cell centered on the origin - prevents trouble
    # if atomic positions are not centered around the origin but displaced
    # far away (remember that translational invariance allows this!)
    #
    ds[:] = RecVecs' * dtau
    ds[:] = ds[:] .- round.(ds)
    dtau0[:] = LatVecs * ds

    # these are estimates of the maximum values of needed integer indices
    nm1 = floor(Int64, norm(bg[:,1]) * rmax ) + 2
    nm2 = floor(Int64, norm(bg[:,2]) * rmax ) + 2
    nm3 = floor(Int64, norm(bg[:,3]) * rmax ) + 2

    for i in -nm1:nm1, j in -nm2:nm2, k in -nm3:nm3
        tt = 0.0
        for ipol in 1:3
            t[ipol] = i*LatVecs[ipol,1] + j*LatVecs[ipol,2] + k*LatVecs[ipol,3] - dtau0[ipol]
            tt = tt + t[ipol]^2
        end
        if (tt <= rmax^2) && (abs(tt) > 1.e-10)
            nrm = nrm + 1
            if nrm > mxr
                error("too many r-vectors, nrm =$(nrm)")
            end
            for ipol in 1:3
                r[ipol,nrm] = t[ipol]
            end
            r2[nrm] = tt
        end
    end
  
    # reorder the vectors in order of increasing magnitude

    if nrm > 1
        @views irr = sortperm(r2[1:nrm])
        @views r2[1:nrm] = r2[irr]
        @views r2[nrm+1:end] .= 0.0
        @views r[:,1:nrm] = r[1:3,irr]
        @views r[:,nrm+1:end] .= 0.0
    end
    return nrm
end

