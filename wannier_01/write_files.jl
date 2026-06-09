
function write_files(NN, LL, nwannier, nband, mu, sigma)
    N1, N2, N3 = NN
    L1, L2, L3 = LL
    #
    #nband = (2L1+1)*(2L2+1)*(2L3+1) #size of the Hamiltonian
    Ks = zeros(Int,3,nband) #list of K vectors
    Ntot = N1*N2*N3
    Ns = [N1,N2,N3]
    R = zeros(ComplexF64,nband,3)
    Rfrac = zeros(ComplexF64,nband,3)
    ind = 1
    for l1 in -L1:L1, l2 in -L2:L2, l3 in -L3:L3
        Ks[:,ind] = [l1, l2, l3]
        R[ind,:]  = 2*pi*[l1/(2*L1+1), l2/(2*L2+1), l3/(2*L3+1)]
        Rfrac[ind,:] = [(l1+L1)/(2*L1+1), (l2+L2)/(2*L2+1), (l3+L3)/(2*L3+1)]
        ind += 1
    end

    win = open("TEMP_free.win","w")
    write(win,"""
    num_bands = $nband
    num_wann = $nwannier

    begin unit_cell_cart
    bohr
    $(2pi) 0 0
    0 $(2pi) 0
    0 0 $(2pi)
    end unit_cell_cart

    mp_grid = $N1 $N2 $N3

    begin kpoints
    """)
    kind = 1
    kind_to_ijk = zeros(Int64,Ntot,3) #kind runs from 1:Ntot, i,j,k from 1:N1,N2,N3
    ijk_to_kind = zeros(Int64,N1,N2,N3)
    kpt_order = [1,2,3] #potentially allow different ordering
    for i in 0:N1-1, j in 0:N2-1, k in 0:N3-1
        write(win, "$(i/N1) $(j/N2) $(k/N3)\n")
        kind += 1
        ijk = [i,j,k]
        big_offset = Ns[kpt_order[1]]*Ns[kpt_order[2]]
        med_offset = Ns[kpt_order[2]]
        kind = ijk[kpt_order[1]]*big_offset + ijk[kpt_order[2]]*med_offset + ijk[kpt_order[3]]+1
        ijk_to_kind[i+1,j+1,k+1] = kind
        kind_to_ijk[kind,:] = [i+1 j+1 k+1]
    end
    write(win, "end kpoints")
    close(win)
    println("Finished writing TEMP_free.win file")

    neighbors = [
        1 0 0;
        0 1 0;
        0 0 1;
       -1 0 0;
        0 -1 0;
        0 0 -1]

    eig = open("TEMP_free.eig","w")
    eig_allk = zeros(ComplexF64,nband,N1*N2*N3)
    ind = 1
    for i in 0:N1-1, j in 0:N2-1, k in 0:N3-1
        ijk = [i,j,k]
        kind = ijk_to_kind[i+1,j+1,k+1]
        eig_k = sort( get_eigvals(ijk./Ns, Ks, NN) )
        eig_allk[:,ind] = eig_k
        ind += 1
        for iband in 1:nband
            write(eig,"$iband $kind $(eig_k[iband])\n")
        end
    end
    close(eig)
    println("Finished writing TEMP_free.eig")

    mmn = open("TEMP_free.mmn","w")
    write(mmn, "Created by free.jl $L1 $L2 $L3\n")
    nneighbors = count(x -> x != 1, Ns)*2
    write(mmn, "$nband $Ntot $nneighbors\n")
    for i in 0:N1-1, j in 0:N2-1, k in 0:N3-1
        ijk = [i,j,k]
        for ib = 1:6
            b = neighbors[ib,:]
            ijkpb_recen = mod.(ijk+b,[N1; N2; N3]) #bring back to BZ
            kind = ijk_to_kind[(ijk.+1)...]
            kindb = ijk_to_kind[(ijkpb_recen.+1)...]
            if kind == kindb
                continue # ignore singleton dimensions
            end
            displ_vec = div.((ijk+b) - ijkpb_recen,Ns)
            write(mmn,"$kind $kindb $(displ_vec[1]) $(displ_vec[2]) $(displ_vec[3])\n")
            ovl = get_ovl(ijk./Ns, (ijkpb_recen)./Ns,displ_vec, Ks, NN, nband)
            for n in 1:nband, m in 1:nband
                write(mmn, string(ovl[m,n]), " 0\n")
            end
        end
    end
    close(mmn)
    println("Finished writing TEMP_free.mmn")

    # Unk and SCDM stuff
    # same for all kpoints
    Amn = zeros(ComplexF64, N1, N2, N3, nband, nwannier)

    gamma = [0.0,0.0,0.0]
    # λ is computed in reduced coordinates
    eig_k = get_eigvals( k_real_to_red(gamma, NN), Ks, NN)
    Unk = exp.(1im*R*Ks)
    occ = (1/2)*erfc.((real(eig_k).-mu)./sigma)

    Qtemp, Rtemp, piv = qr(diagm(0=>occ)*Unk', Val(true))

    cols = piv[1:nwannier]
    min_S = Inf
    for i in 0:N1-1, j in 0:N2-1, k in 0:N3-1
        kpt = [i,j,k]./Ns

        phase = exp.( 1im*R[cols,:]*k_red_to_real(kpt, NN) )

        eig_k = get_eigvals(kpt, Ks, NN)
        n_to_K = sortperm(eig_k)
        K_to_n = invperm(n_to_K)
        eig_k = eig_k[n_to_K]
        # Note that the order of Ks is kpt dependent
        Unk = zeros(ComplexF64,nband,nband)
        for m = 1 : nband
            Km = n_to_K[m]
            Unk[:,m] = exp.(1im*R*(Ks[:,Km]))
        end
        occ = (1/2)*erfc.((real(eig_k).-mu)./sigma)

        Y = diagm(0=>phase)*Unk[cols,:]*diagm(0=>occ)
        Y = Y'

        U,S,V = svd(Y)
        Amn[i+1,j+1,k+1,:,:] = U*V'
        min_S = min(min_S, minimum(S))
    end
    println("Min singular value: $min_S")

    for file in ("TEMP_free.amn", "TEMP_free.SCDM.amn")
        out = open(file,"w")
        write(out, "Created by free.jl ", string(Dates.now()), "\n")
        write(out, "$(nband) $(N1*N2*N3) $(nwannier)\n")
        for i in 1:N1, j in 1:N2, k in 1:N3
            kind = ijk_to_kind[i,j,k]
            for n in 1:nwannier, m in 1:nband
                coeff = Amn[i,j,k,m,n]
                write(out, "$m $n $kind $(real(coeff)) $(imag(coeff))\n")
            end
        end
        close(out)
    end
    return
end
