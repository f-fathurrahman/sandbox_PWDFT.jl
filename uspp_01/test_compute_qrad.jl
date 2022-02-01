using PWDFT

function create_atoms_N2H4()
    atoms = Atoms(xyz_string="""
    6

    N       5.94821400       6.81171100       5.22639100
    N       5.94821400       5.37379300       5.22639100
    H       6.15929600       7.18550400       6.15196500
    H       5.00000000       7.09777800       5.00000000
    H       5.73713200       5.00000000       6.15196500
    H       6.89642800       5.08772600       5.00000000
    """, LatVecs=gen_lattice_sc(16.0))
    return atoms
end

function compute_qrad(atoms, pw, pspots)
    CellVolume = pw.CellVolume
    Nspecies = atoms.Nspecies
    prefr = 4π/CellVolume
    ecutrho = pw.ecutrho

    ndm = 0
    for isp in 1:Nspecies
        println("kkbeta = ", pspots[isp].kkbeta)
        if ndm < pspots[isp].kkbeta
            ndm = pspots[isp].kkbeta
        end
    end
    println("ndm = ", ndm)

    qnorm = 0.0 # XXX HARDCODED, no k-points norm of (q + k) ?
    dq = 0.01 # XXX HARDCODED
    cell_factor = 1.0 # hardcoded

    #ndm = max( upf(:)%kkbeta )
    #nqxq = INT( ( (SQRT(ecutrho) + qnorm) / dq + 4) * cell_factor )
    nqxq = round(Int64, sqrt(2*ecutrho)/dq + 4) # convert to Ry
    println("nqxq = ", nqxq)
end




function main()
    atoms = create_atoms_N2H4()
    println(atoms)

    pw = PWGrid(20.0, atoms.LatVecs, dual=5.0)
    println(pw)
    
    Nspecies = atoms.Nspecies
    pspfiles = [
        "/home/efefer/pseudo/GBRV_LDA/n_lda_v1.2.uspp.F.UPF2",
        "/home/efefer/pseudo/GBRV_LDA/h_lda_v1.4.uspp.F.UPF2"
    ]
    pspots = Array{PsPot_UPF}(undef,Nspecies)
    for isp = 1:Nspecies
        pspots[isp] = PsPot_UPF( pspfiles[isp] )
        PWDFT._build_prj_interp_table!( pspots[isp], pw )
    end

    compute_qrad(atoms, pw, pspots)

end

main()