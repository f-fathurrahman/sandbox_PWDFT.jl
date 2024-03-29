function calc_forces_Ps_loc( Ham )
    F_Ps_loc = zeros(Float64, 3, Ham.atoms.Natoms)
    calc_forces_Ps_loc!( Ham, F_Ps_loc )
    return F_Ps_loc
end

function calc_forces_Ps_loc!( Ham::Hamiltonian, F_Ps_loc::Array{Float64,2} )
    calc_forces_Ps_loc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe, F_Ps_loc )
    return
end


function calc_forces_Ps_loc!(
    atoms::Atoms, pw::PWGrid, pspots::Array{PsPot_GTH,1},
    Rhoe::Array{Float64,2},
    F_Ps_loc
)
    
    Ng = pw.gvec.Ng
    Natoms = atoms.Natoms

    G = pw.gvec.G
    G2 = pw.gvec.G2

    idx_g2r = pw.gvec.idx_g2r
    
    Npoints = prod(pw.Ns)
    
    Rhoe_tot = zeros(Float64,Npoints)
    Nspin = size(Rhoe,2)
    for ispin in 1:Nspin, ip in 1:Npoints
        Rhoe_tot[ip] = Rhoe_tot[ip] + Rhoe[ip,ispin]
    end

    RhoeG = R_to_G(pw, Rhoe_tot)/Npoints  # this normalization is required

    Ω = pw.CellVolume

    atpos = atoms.positions
    
    for ia in 1:Natoms
        isp = atoms.atm2species[ia]
        psp = pspots[isp]
        # G=0 contribution is zero, G[:,ig=0] = [0, 0, 0]
        for ig in 2:Ng
            @views GX = dot(atpos[:,ia], G[:,ig])
            Sf = cos(GX) - im*sin(GX)
            VlocG = eval_Vloc_G( psp, G2[ig] )
            ip = idx_g2r[ig]
            @views F_Ps_loc[:,ia] .+= real( im*G[:,ig] * VlocG * conj(RhoeG[ip])*Sf )
        end
    end

    return
end



function calc_forces_Ps_loc_finite_diff(
    atoms::Atoms, pw::PWGrid, pspots::Array{PsPot_GTH,1},
    Rhoe::Array{Float64,2}
)
    Δ = 0.001

    pos_orig = copy(atoms.positions)
    Natoms = atoms.Natoms

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    Rhoe_tot = zeros(Float64,Npoints)
    Nspin = size(Rhoe)[2]
    for ispin = 1:Nspin
        for ip = 1:Npoints
            Rhoe_tot[ip] = Rhoe_tot[ip] + Rhoe[ip,ispin]
        end
    end
    @printf("Integ Rhoe_tot = %18.10f\n", sum(Rhoe_tot)*pw.CellVolume/prod(pw.Ns))

    F_Ps_loc = zeros(3,Natoms)
    V_Ps_loc = zeros(Float64, Npoints)

    for ia = 1:Natoms
    for idir = 1:3
        
        # set to original positions
        atoms.positions[:,:] = pos_orig[:,:]
        
        atoms.positions[idir,ia] = pos_orig[idir,ia] + 0.5*Δ        
        V_Ps_loc = calc_V_Ps_loc(atoms, pw, pspots)
        Eplus = dot( V_Ps_loc, Rhoe_tot ) * dVol
        
        atoms.positions[idir,ia] = pos_orig[idir,ia] - 0.5*Δ
        V_Ps_loc = calc_V_Ps_loc(atoms, pw, pspots)
        Eminus = dot( V_Ps_loc, Rhoe_tot ) * dVol

        F_Ps_loc[idir,ia] = -(Eplus - Eminus)/Δ
    end
    end

    atoms.positions[:,:] = pos_orig[:,:]

    return F_Ps_loc
end


function calc_V_Ps_loc(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Array{PsPot_GTH,1}
)

    Npoints = prod(pw.Ns)
    Nspecies = atoms.Nspecies    
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    CellVolume = pw.CellVolume
    idx_g2r = pw.gvec.idx_g2r

    Vg = zeros(ComplexF64,Npoints)
    V_Ps_loc = zeros(Float64,Npoints)

    strf = calc_strfact( atoms, pw )

    for isp = 1:Nspecies
        psp = pspots[isp]
        for ig = 1:Ng
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] ) / CellVolume
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    return V_Ps_loc
end