# Update the Ham.atoms.positions, local and nonlocal pseudopotentials
function update_positions!( Ham::Hamiltonian, dr::Array{Float64,2} )
    Ham.atoms.positions[:] = Ham.atoms.positions[:] + dr[:]

    # Make sure that the symmetry is not changed
    # Probably this is not mandatory, but the some fields
    # of Hamiltonian need to be changed.
    println("Nsyms = ", Ham.sym_info.Nsyms)
    if Ham.sym_info.Nsyms > 1
        println("Nsyms = ", Ham.sym_info.Nsyms)
        new_sym_info = SymmetryInfo(Ham.atoms)
        @assert new_sym_info.Nsyms == Ham.sym_info.Nsyms
    end

    # pwgrid and kpoints should not change.

    strf = calc_strfact( Ham.atoms, Ham.pw )

    pspots = Ham.pspots
    G2_shells = Ham.pw.gvec.G2_shells
    idx_g2shells = Ham.pw.gvec.idx_g2shells
    Npoints = prod(Ham.pw.Ns)
    CellVolume = Ham.pw.CellVolume
    Ng = Ham.pw.gvec.Ng
    Nspecies = Ham.atoms.Nspecies
    idx_g2r = Ham.pw.gvec.idx_g2r

    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)

    for isp in 1:Nspecies
        psp = pspots[isp]
        # FIXME: This is not really depends on ionic positions
        # Probably we should save Vgl in Hamiltonian
        for igl in 1:Ngl
            Vgl[igl] = eval_Vloc_G( psp, G2_shells[igl] )
        end
        #
        for ig in 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            Vg[ip] = strf[ig,isp] * Vgl[igl]
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(Ham.pw, Vg) ) * Npoints / CellVolume
    end
    Ham.potentials.Ps_loc[:] = V_Ps_loc[:]

    Ham.pspotNL = PsPotNL( Ham.atoms, Ham.pw, pspots, check_norm=false )

    return
end