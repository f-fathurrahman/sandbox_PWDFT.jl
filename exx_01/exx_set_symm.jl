function exx_set_symm( sym_info, Ns, Nsx )

    Nsyms = sym_info.Nsyms

    nr1, nr2, nr3 = Ns
    nr1x, nr2x, nr3x = Nsx  

    nxxs = nr1x*nr2x*nr3x
    rir = zeros(Int64, nxxs, Nsyms)
    
    ftau = zeros(Int64, 3, Nsyms)
    s_scaled = zeros(Int64, 3, 3, Nsyms)

    scale_sym_ops!(sym_info, nr1, nr2, nr3, s_scaled, ftau)
    for isym in 1:Nsyms
        for k in 1:nr3, j in 1:nr2, i in 1:nr1
            ri, rj, rk = rotate_grid_point(
                s_scaled[:,:,isym], ftau[:,isym],
                i, j, k, nr1, nr2, nr3 
            )
            ir = i + (j-1)*nr1x + (k-1)*nr1x*nr2x
            rir[ir,isym] = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
        end
    end
    return rir
end