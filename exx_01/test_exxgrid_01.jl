using Infiltrator
using Printf
using LinearAlgebra: norm
using PWDFT

# will modify v
function cryst_to_cart!(nv, v, trmat, iflag)
    @assert iflag in [-1,1]
    vau = zeros(Float64, 3)
    for nv in 1:nv
        if iflag == 1
            for k in 1:3
                vau[k] = trmat[k,1]*v[1,nv] + trmat[k,2]*v[2,nv] + trmat[k,3]*v[3,nv]
            end
        else
            for k in 1:3
                vau[k] = trmat[1,k]*v[1,nv] + trmat[2,k]*v[2,nv] + trmat[3,k]*v[3,nv]
            end
        end
        for k in 1:3
            v[k,nv] = vau[k]
        end
    end
    return
end


function exx_qgrid_init!(
    pw, nqs, nq1, nq2, nq3,
    max_nk, temp_nkqs, temp_xkq, temp_index_ikq, dxk
)
    SMALL_Q = 1e-6 #XXX FIXME: Hardcoded

    Nkpt = pw.gvecw.kpoints.Nkpt

    new_ikq = zeros(Int64, max_nk)
    index_xkq = zeros(Int64, Nkpt, nqs)
    # nqs is computed with current nq1, nq2, nq3

    # should be returned
    nkqs = 0 

    sxk = zeros(Float64, 3)
    xk_cryst = zeros(Float64, 3)

    # define the q-mesh step-sizes
    dq1 = 1.0/nq1
    dq2 = 1.0/nq2
    dq3 = 1.0/nq3
    for ik in 1:Nkpt
        # go to crystalline coordinates
        xk_cryst[:] = pw.gvecw.kpoints.k[:,ik]
        cryst_to_cart!(1, xk_cryst, pw.LatVecs/(2π), -1)
        #
        iq = 0
        for iq1 in 1:nq1, iq2 in 1:nq2, iq3 in 1:nq3
            sxk[1] = xk_cryst[1] + (iq1-1)*dq1
            sxk[2] = xk_cryst[2] + (iq2-1)*dq2
            sxk[3] = xk_cryst[3] + (iq3-1)*dq3
            iq += 1
            xk_not_found = true
            for ikq in 1:temp_nkqs
                if xk_not_found
                    dxk[:] = sxk[:] - temp_xkq[:,ikq] - round.(Int64, sxk[:] - temp_xkq[:,ikq])
                    if all( abs.(dxk) .<= SMALL_Q )
                        xk_not_found = false
                        if new_ikq[ikq] == 0
                            nkqs += 1
                            temp_index_ikq[nkqs] = ikq
                            new_ikq[ikq] = nkqs
                        end
                        index_xkq[ik,iq] = new_ikq[ikq]
                    end # if
                end # if
            end # for ikq
            #
            if xk_not_found
                return nkqs, index_xkq
            end
        end
    end #
    return nkqs, index_xkq
end



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

    for isym in 1:Nsyms
        println("\nisym = ", isym)
        for ik in 1:Nkpt
            # go to crystalline coordinates
            xk_cryst[:] = pw.gvecw.kpoints.k[:,ik]
            #println()
            #@printf("%4d in bg   = [%18.10f %18.10f %18.10f]\n", ik, alat*xk_cryst...)
            cryst_to_cart!( 1, xk_cryst, pw.LatVecs/(2π), -1)
            #xk_cryst[:] = pw.LatVecs * xk_cryst / (2π) same as this
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

    @infiltrate

    return

end

