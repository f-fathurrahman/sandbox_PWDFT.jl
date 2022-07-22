function update_from_rhoe!(Ham, Rhoe, RhoeG)

    Ham.rhoe[:,:] = Rhoe[:,:] # Need copy?

    # FIXME: also set Ham.potentials.XC and Ham.potentials.Hartree
    Exc, Evtxc = _add_V_xc!( Ham, Rhoe, RhoeG )
    Ehartree = _add_V_Hartree!( Ham, Rhoe, RhoeG )
    
    Veff = Ham.potentials.Total
    println("After add_V_Hartree! sum abs Veff = ", sum(Veff))
    println("sum potentials.Ps_loc = ", sum(Ham.potentials.Ps_loc))

    Nspin = Ham.electrons.Nspin
    for ispin in 1:Nspin
        Ham.potentials.Total[:,ispin] .+= Ham.potentials.Ps_loc[:]
    end

    if Ham.pw.using_dual_grid
        dense_to_smooth!( Ham.pw, Ham.potentials.Total, Ham.potentials.TotalSmooth )
    end

    return # energies?
end



function _add_V_xc!(Ham, Rhoe, RhoeG)

    Nspin = Ham.electrons.Nspin
    
    @assert Nspin == 1

    Npoints = size(Rhoe, 1)
    epsxc = zeros(Float64, Npoints)
    Vxc = Ham.potentials.XC

    # XC potential
    # VWN is the default
    if Ham.rhoe_core == nothing
        epsxc[:], Vxc[:,1] = calc_epsxc_Vxc_VWN( Ham.xc_calc, Rhoe[:,1] )
    else
        epsxc[:], Vxc[:,1] = calc_epsxc_Vxc_VWN( Ham.xc_calc, Rhoe[:,1] + Ham.rhoe_core )
    end

    dVol = Ham.pw.CellVolume / Npoints
    # Calculate etxc and vtxc
    Exc = sum(epsxc .* (Rhoe[:,1] + Ham.rhoe_core))*dVol
    Evtxc = sum(Vxc[:,1] .* Rhoe)*dVol # Evtxc does not include rhoe_core

    Ham.potentials.Total[:,1] += Vxc[:,1] # Update

    return Exc, Evtxc
end


# Note that RhoeG is already in FFT grid
function _add_V_Hartree!(Ham, Rhoe, RhoeG)

    pw = Ham.pw
    gvec = pw.gvec

    G2 = gvec.G2
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    Npoints = prod(pw.Ns) # dense

    VhG = zeros(ComplexF64, Npoints)
    Ehartree = 0.0
    # skip ig = 1, it is set to zero
    for ig in 2:Ng
        ip = idx_g2r[ig]
        Ehartree = Ehartree + 2π*( real(RhoeG[ip])^2 + imag(RhoeG[ip])^2 )/G2[ig]
        VhG[ip] = 4π * RhoeG[ip]/G2[ig]
    end

    Ehartree *= Ham.pw.CellVolume

    G_to_R!(pw, VhG)
    VhG[:] *= Npoints # XXX: scale by Npoints
    Ham.potentials.Hartree[:] = real(VhG) # update
    Ham.potentials.Total[:,1] += real(VhG)

    return Ehartree
end