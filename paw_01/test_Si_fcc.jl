using PWDFT
using LinearAlgebra

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials")

function create_atoms_Si_fcc()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    return atoms
end

function create_Ham_Si_fcc_paw_pslib()
    atoms = create_atoms_Si_fcc()

    # Initialize Hamiltonian
    pspots = [
        PsPot_UPF(joinpath(DIR_PWDFT, "pseudopotentials",
            "PSLIB_US_PAW_LDA", "Si.pz-n-kjpaw_psl.1.0.0.UPF"))
    ]
    # FIXME: PSLIB_US_PAW_LDA is not included in the repo
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    options.meshk = [3,3,3]
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end


function PAW_Ehxc_energy!( Ham::Hamiltonian )
    return PAW_Ehxc_energy!(
        Ham.atoms, Ham.pspots, Ham.pspotNL, Ham.xc_calc,
        Ham.pspotNL.becsum,
        Ham.pspotNL.paw.E_paw_cmp
    )
end


# Similar to PAW_potential, but without ddd_paw calculation
function PAW_Ehxc_energy!(
    atoms::Atoms,
    pspots::Vector{PsPot_UPF},
    pspotNL::PsPotNL_UPF,
    xc_calc,
    becsum, e_cmp
)

    fill!(e_cmp, 0.0)
    energy_tot = 0.0

    Nspin = size(becsum, 3)

    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nh = pspotNL.nh

    # Preallocate work arrays (for each species)
    # XXX: Use Nrmesh max ?
    v_lm_s = Vector{Array{Float64,3}}(undef,atoms.Nspecies)
    savedv_lm_s = Vector{Array{Float64,3}}(undef,atoms.Nspecies)
    rho_lm_s = Vector{Array{Float64,3}}(undef,atoms.Nspecies)
    for isp in 1:Nspecies
        Nrmesh = pspots[isp].Nr
        l2 = (pspots[isp].lmax_rho + 1)^2
        v_lm_s[isp] = zeros(Float64, Nrmesh, l2, Nspin)
        rho_lm_s[isp] = zeros(Float64, Nrmesh, l2, Nspin)
    end


    # Begin loop over all atoms
    for ia in 1:Natoms

        isp = atm2species[ia]
        # Skip this atoms if it is not using PAW
        if !pspots[isp].is_paw
            continue
        end

        Nrmesh = pspots[isp].Nr
        l2 = (pspots[isp].lmax_rho + 1)^2
        kkbeta = pspots[isp].kkbeta

        # We might need to call GC manually to free these arrays
        v_lm = v_lm_s[isp]
        rho_lm = rho_lm_s[isp]

        # All-electron contribution
        for AE in [true, false]
            if AE
                i_what = 1
                sgn = 1
            else
                i_what = 2
                sgn = -1
            end
            # sgn: sign for energy summation

            PAW_rho_lm!(AE, ia, atoms, pspots, pspotNL, becsum, rho_lm)

            # Hartree term
            @views energy = PAW_h_potential!( ia, atoms, pspots, rho_lm, v_lm[:,:,1] )
            energy_tot += sgn*energy
            e_cmp[ia,1,i_what] = sgn*energy # Hartree, all-electron

            # XC term
            @views energy = PAW_xc_potential!( AE, ia, atoms, pspots, pspotNL, xc_calc, rho_lm, v_lm )
            energy_tot += sgn*energy
            e_cmp[ia,2,i_what] = sgn*energy # XC, all-electron

        end # AE, PS
    end # loop over all atoms

    return energy_tot
end


function main()
    Ham = create_Ham_Si_fcc_paw_pslib()
    psiks = rand_BlochWavefunc(Ham)
    
    use_smearing = false
    electrons_scf!(Ham, psiks)

    energies = calc_energies(Ham, psiks)
    Ehxc = PAW_Ehxc_energy!(Ham)
    if use_smearing
        energies.mTS = mTS
    end
    println("\nUsing original formula for total energy")
    println(energies, use_smearing=use_smearing)
    println("Ehxc = ", Ehxc)
    println("Total energy (in Ha) = ", sum(energies) + Ehxc)
    println("! Total energy (in Ry) = ", 2*(sum(energies) + Ehxc))

end

main()
