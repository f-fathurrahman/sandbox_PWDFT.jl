# Update the Ham.atoms.positions, local and nonlocal pseudopotentials
function update_positions!( Ham::Hamiltonian, new_pos::Array{Float64,2} )
    Ham.atoms.positions[:] = new_pos[:]

    println("New poisitons: ")
    display(Ham.atoms.positions'); println();

    # pwgrid and kpoints should not change.

    atoms = Ham.atoms
    pw = Ham.pw
    strf = calc_strfact( atoms.positions, atoms.Nspecies, atoms.atm2species, pw.gvec.G )

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
        #
        for igl = 1:Ngl
            Vgl[igl] = eval_Vloc_G( psp, G2_shells[igl] )
        end
        #
        Vg[1] = strf[1,isp]*Vgl[1] # G=(0,0,0)
        for ig in 2:Ng
            #
            igl = idx_g2shells[ig]
            #
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * Vgl[igl]
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints / CellVolume
    end
    Ham.potentials.Ps_loc[:] = V_Ps_loc[:]

    Ham.pspotNL = PsPotNL( Ham.atoms, Ham.pw, pspots )

    return
end