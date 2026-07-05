using Infiltrator
using Printf
using LinearAlgebra: norm, inv
using PWDFT

includet("cryst_to_cart.jl")
includet("exx_grid_check.jl")
includet("exx_qgrid_init.jl")

function debug_main()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nsyms = Ham.sym_info.Nsyms
    nq1 = 1
    nq2 = 1
    nq3 = 1

    max_nk = Nkpt * min(48, 2*Nsyms)
    temp_index_xk = zeros(Int64, max_nk)
    temp_index_sym = zeros(Int64, max_nk)
    temp_index_ikq = zeros(Int64, max_nk)
    temp_xkq = zeros(Float64, 3, max_nk)
    
    pw = Ham.pw

    # find all k-points equivalent by symmetry to the points in the k-list
    temp_nkqs = 0
    xk_cryst = zeros(Float64, 3)
    sxk = zeros(Float64, 3)
    s = Ham.sym_info.s
    SMALL_Q = 1e-6
    noncolin = Ham.electrons.noncollinear
    domag = Ham.electrons.domag
    dxk = zeros(Float64, 3)

    for isym in 1:Nsyms
        println("\nisym = ", isym)
        for ik in 1:Nkpt
            # go to crystalline coordinates
            xk_cryst[:] = pw.gvecw.kpoints.k[:,ik]
            println()
            @printf("%4d in bg   = [%18.10f %18.10f %18.10f]\n", ik, xk_cryst...)
            #cryst_to_cart!( 1, xk_cryst, pw.LatVecs/(2π), -1)
            xk_cryst[:] = inv(pw.RecVecs)*xk_cryst[:] # This is also can be used, but using inverse
            #cryst_to_cart!( 1, xk_cryst, pw.RecVecs, 1) # cannot use RecVecs
            #cryst_to_cart!( 1, xk_cryst, pw.LatVecs, -1) # not working without 2π
            #xk_cryst[:] = pw.LatVecs * xk_cryst / (2π) same as this?
            #println("$ik $xk_cryst")
            @printf("%4d in cart = [%18.10f %18.10f %18.10f]\n", ik, xk_cryst...)
            # rotate with this sym.op.
            sxk[:] = s[:,1,isym]*xk_cryst[1] + s[:,2,isym]*xk_cryst[2] + s[:,3,isym]*xk_cryst[3]
            # add sxk to the auxiliary list IF it is not already present
            xk_not_found = true
            # do-loop skipped the first time because temp_nkqs == 0
            for ikq in 1:temp_nkqs
                if xk_not_found
                    dxk[:] = sxk[:] - temp_xkq[:,ikq] - round.(Int64, sxk[:] - temp_xkq[:,ikq])
                    if all( abs.(dxk) .<= SMALL_Q )
                        xk_not_found = false
                    end
                end
            end
            #
            if xk_not_found
                temp_nkqs += 1
                temp_xkq[:,temp_nkqs] = sxk[:]
                temp_index_xk[temp_nkqs] = ik
                temp_index_sym[temp_nkqs] = isym
            end
            #
            sxk[:] = -sxk[:]
            xk_not_found = true
            for ikq in 1:temp_nkqs
                if xk_not_found
                    dxk[:] = sxk[:] - temp_xkq[:,ikq] - round.(Int64, sxk[:] - temp_xkq[:,ikq])
                    if all( abs.(dxk) .<= SMALL_Q )
                        xk_not_found = false
                    end
                end
            end
            if xk_not_found && !(noncolin && domag)
                temp_nkqs += 1
                temp_xkq[:,temp_nkqs] = sxk[:]
                temp_index_xk[temp_nkqs] = ik
                temp_index_sym[temp_nkqs] = -isym
            end
        end
    end
    println("temp_nkqs = ", temp_nkqs)
    
    nqs = nq1*nq2*nq3
    nkqs, index_xkq = exx_qgrid_init!(
        pw, nqs, nq1, nq2, nq3,
        max_nk, temp_nkqs, temp_xkq, temp_index_ikq, dxk
    )

    nspin_lsda = 1 # XXX hardcoded

    xkq_collect = zeros(Float64, 3, nspin_lsda*nkqs)
    index_xk = zeros(Int64, nspin_lsda*nkqs)
    index_sym = zeros(Int64, nspin_lsda*nkqs)
    for ik in 1:nkqs
        ikq = temp_index_ikq[ik]
        xkq_collect[:,ik] = temp_xkq[:,ikq]
        index_xk[ik] = temp_index_xk[ikq]
        index_sym[ik] = temp_index_sym[ikq]
    end

    #println("xkq_collect in fractional:")
    #for ik in 1:nkqs
    #    @printf("[%8.5f %8.5f %8.5f] %4d %4d\n",
    #        xkq_collect[1,ik], xkq_collect[2,ik], xkq_collect[3,ik], 
    #        index_xk[ik], index_sym[ik])
    #end
    xkq_collect[:,:] = pw.RecVecs*xkq_collect[:,:]
    println("xkq_collect in Cartesian:")
    for ik in 1:nkqs
        @printf("[%8.5f %8.5f %8.5f] %4d %4d\n",
            xkq_collect[1,ik], xkq_collect[2,ik], xkq_collect[3,ik], 
            index_xk[ik], index_sym[ik])
    end
    #xkq_collect[:,:] = inv(pw.RecVecs)*xkq_collect[:,:]
    #println("xkq_collect in fractional (again):")
    #for ik in 1:nkqs
    #    @printf("[%8.5f %8.5f %8.5f] %4d %4d\n",
    #        xkq_collect[1,ik], xkq_collect[2,ik], xkq_collect[3,ik], 
    #        index_xk[ik], index_sym[ik])
    #end

    exx_grid_check(
        pw, Ham.sym_info,
        nq1, nq2, nq3, index_xk, index_xkq, index_sym
    )

    qnrm = 0.0
    for iq in 1:nkqs, ik in 1:Nkpt
        qnrm = max(qnrm, sqrt( sum( (pw.gvecw.kpoints.k[:,ik] - xkq_collect[:,iq]).^2) ))
    end
    println("qnrm = ", qnrm)

    @infiltrate

    return

end

