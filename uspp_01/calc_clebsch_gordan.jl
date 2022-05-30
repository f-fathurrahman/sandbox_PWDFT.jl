# From aainit of QE

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

# lli = lmaxkb + 1
function calc_clebsch_gordan( lmaxkb::Int64 )

    lli = lmaxkb + 1

    # Those are parameters (HARDCODED)
    lmaxx  = 3         # max non local angular momentum (l=0 to lmaxx)      
    lqmax = 2*lmaxx + 1  # max number of angular momenta of Q

    # maximum number of combined angular momentum
    nlx = (lmaxx + 1)^2
    # maximum magnetic angular momentum of Q
    mx = 2*lqmax - 1

    llx = (2*lli - 1)^2

    @assert (2*lli-1) <= lqmax
    @assert lli >= 0
    @assert (lli*lli) <= nlx

    println("lli = ", lli)
    println("llx = ", llx)

    r = zeros(Float64, 3, llx)
    for ir in 1:llx
        c = 2.0*rand() - 1.0
        s = sqrt(1.0 - c^2)
        ϕ = 2π * rand()
        r[1,ir] = s * cos(ϕ)
        r[2,ir] = s * sin(ϕ)
        r[3,ir] = c
    end

    Ylm = zeros(Float64, llx, llx)
    # generate the real spherical harmonics for the array: ylm(ir,lm)
    _lmax = round(Int64, sqrt(llx) - 1)
    Ylm_real_qe!(_lmax, r, Ylm)
    #for i in 1:llx
    #    @views Ylm_real_qe!(_lmax, r[:,i], Ylm[i,:])
    #end

    Ylminv = inv(Ylm)

    # Clebsch-Gordan coefficients for spherical harmonics
    ap = zeros(Float64, lqmax*lqmax, nlx, nlx)
    # for each pair of combined momenta lm(1),lm(2): 
    lpx = zeros(Int64, nlx, nlx)      # maximum combined angular momentum LM
    lpl = zeros(Int64, nlx, nlx, mx)  # list of combined angular momenta  LM

    # for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
    println()
    for li in 1:lli*lli
        println()
        for lj in 1:lli*lli
            lpx[li,lj] = 0
            for l in 1:llx
                ap[l,li,lj] = _compute_ap(l, li, lj, llx, Ylm, Ylminv)
                if abs(ap[l,li,lj]) > 1e-3
                    lpx[li,lj] = lpx[li,lj] + 1   # increment
                    if lpx[li,lj] > mx
                        println("Error: mx dimension too small: lpx(li,lj)")
                        error("Error in calc_clebsch_gordan")
                    end
                    lpl[li, lj, lpx[li,lj]] = l
                    @printf("%4d %4d %4d %18.10f\n", l, li, lj, ap[l,li,lj])
                end # if
            end
        end
    end

    return ap

end



function _compute_ap(l,li,lj,llx,ylm,mly)
    res = 0.0
    for i in 1:llx
        res += mly[l,i]*ylm[i,li]*ylm[i,lj]
    end
    return res
end
