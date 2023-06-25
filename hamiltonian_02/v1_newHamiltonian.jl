function newHamiltonian(
    atoms::Atoms,
    pspots::Vector{Tpsp},
    ecutwfc::Float64,
    options::HamiltonianOptions
) where Tpsp

    if options.use_symmetry == false
        sym_info = SymmetryInfo()
    else
        sym_info = SymmetryInfo(atoms)
    end

    @assert options.dual >= 4.0

    # kpoints
    if options.kpoints == nothing
        if options.kpts_str == ""
            # automatic generation of kpoints
            kpoints = KPoints( atoms,
                options.meshk, options.shiftk,
                sym_info.s, time_reversal=options.time_reversal )
        else
            # use the given kpoints
            kpoints = kpoints_from_string( atoms, options.kpts_str )
        end
    else
        @assert typeof(kpoints) == KPoints
    end

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc,
        atoms.LatVecs, kpoints=kpoints,
        Ns_=options.Ns,
        dual=options.dual
    )

    Nspecies = atoms.Nspecies
    @assert size(pspots,1) == Nspecies

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    idx_g2shells = pw.gvec.idx_g2shells
    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)


    #
    # Initialize pseudopotentials and local ionic potentials
    #
    is_psp_using_nlcc = zeros(Bool,Nspecies) # by default we don't use any NLCC    
    V_of_0 = 0.0
    #
    for isp in 1:Nspecies
        psp = pspots[isp] # shortcut
        if psp.is_nlcc
            is_psp_using_nlcc[isp] = true
        end
        eval_Vloc_G!( psp, G2_shells, Vgl )
        for ig in 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            Vg[ip] += strf[ig,isp] * Vgl[igl] / CellVolume
        end
    end
    V_of_0 += real(Vg[1]) # using PWSCF convention
    V_Ps_loc = real(G_to_R(pw, Vg))
    V_Ps_loc *= Npoints # Rescale using PWSCF convention

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_xc = zeros( Float64, Npoints, options.Nspin )
    V_loc_tot = zeros( Float64, Npoints, options.Nspin )
    if pw.using_dual_grid
        # We initialize smooth local potential here (total)
        potentials = Potentials(
            V_Ps_loc, V_Hartree, V_xc, V_loc_tot,
            zeros(Float64, prod(pw.Nss), options.Nspin),
            zeros(Float64, Npoints, options.Nspin)
        )
    else
        potentials = Potentials(
            V_Ps_loc, V_Hartree, V_xc, V_loc_tot,
            nothing,
            zeros(Float64, Npoints, options.Nspin)
        )
    end


    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, options.Nspin )
    #
    # Initialize core electron density for NLCC if needed
    #
    if any(is_psp_using_nlcc)
        rhoe_core = zeros(Float64, Npoints, 1) # FIXME: Nspin=1
        calc_rhoe_core!(atoms, pw, pspots, rhoe_core)
    else
        rhoe_core = nothing
    end


    if options.extra_states > -1
        electrons = Electrons( atoms, pspots, Nspin=options.Nspin, Nkpt=kpoints.Nkpt,
            Nstates_empty=extra_states )
    elseif options.Nstates > -1
        electrons = Electrons( atoms, pspots, Nspin=options.Nspin, Nkpt=kpoints.Nkpt,
            Nstates=options.Nstates )
    elseif (options.Nstates == -1) && (options.extra_states == -1)
        # Default value for Nstates and Nstates_empty
        # Nstates will be calculated automatically
        electrons = Electrons( atoms, pspots, Nspin=options.Nspin, Nkpt=kpoints.Nkpt )
    else
        error("Error in initializing instance of Electrons")
    end


    # FIXME: Make parametric PsPotNL ?
    # NL pseudopotentials
    are_using_upf = zeros(Bool, Nspecies)
    for isp in 1:Nspecies
        pspfile = pspots[isp].pspfile
        are_using_upf[isp] = PWDFT.is_using_extension_upf(pspfile)
    end
    if all(are_using_upf)
        pspotNL = PsPotNL_UPF(atoms, pw, pspots)
    elseif all(.!are_using_upf)
        pspotNL = PsPotNL( atoms, pw, pspots, check_norm=false )
    else
        error("Not supporting mixed pseudopotential types: GTH and UPF")
    end
    # XXX For the moment we treat GTH and UPF as different.
    # We also do not support mixing UPF and GTH pspots
    # Ideally we should have one PsPotNL type only.


    atoms.Zvals = get_Zvals( pspots )

    ik = 1
    ispin = 1

    if sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( atoms, pw, sym_info )
    else
        rhoe_symmetrizer = RhoeSymmetrizer() # dummy rhoe_symmetrizer
    end

    if options.use_xc_internal
        xc_calc = XCCalculator()
    else
        # Using Libxc is the default
        if options.xcfunc == "SCAN"
            xc_calc = LibxcXCCalculator(is_metagga=true, Npoints=Npoints, Nspin=options.Nspin)
        elseif options.xcfunc == "PBE"
            xc_calc = LibxcXCCalculator(x_id=101, c_id=130)
        else
            xc_calc = LibxcXCCalculator()
        end
    end

    return Hamiltonian( pw, potentials, energies, rhoe, rhoe_core,
                        electrons, atoms, sym_info, rhoe_symmetrizer,
                        pspots, pspotNL, options.xcfunc, xc_calc, ik, ispin )
end
