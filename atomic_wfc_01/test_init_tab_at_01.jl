using Revise

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials")

#=
Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");
=#

includet("init_tab_at.jl")

function test_main(; filename=nothing, do_export_data=false)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    tab_at = init_tab_at(Ham.pspots[1], Ham.pw)
end

function create_Ham_structure_01()
    atoms_tuple = (
        Natoms = 2,
        Nspecies = 1,
        positions = [0.0 2.7079775377810678; 0.0 2.7079775377810678; 0.0 2.7079775377810678],
        atm2species = [1, 1],
        atsymbs = ["Fe", "Fe"],
        SpeciesSymbols = ["Fe"],
        LatVecs = [5.4159550755621355 0.0 0.0; 0.0 5.4159550755621355 0.0; 0.0 0.0 5.4159550755621355],
        Zvals = [16.0],
        masses = [0.0]
    );
    atoms = Atoms(atoms_tuple...);
    
    pspfiles = [ joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "Fe.upf") ];
    ecutwfc = 20.0;
    options_tuple = (
        dual = 4.0, Nspin_wf = 2, Nspin_dens = 2, meshk = [3, 3, 3],
        shiftk = [0, 0, 0],
        time_reversal = true, Ns = (0, 0, 0), kpoints = nothing, kpts_str = nothing,
        xcfunc = "VWN",
        use_xc_internal = false, extra_states = nothing,
        Nstates = 20, use_symmetry = true, use_smearing = true,
        smearing_kT = 0.001, starting_magn = [0.4], angle1 = nothing,
        angle2 = nothing, lspinorb = false, noncollinear = false
    );
    options = HamiltonianOptions(options_tuple...);
    pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
    for isp in 1:atoms.Nspecies
        pspots[isp] = PsPot_UPF(pspfiles[isp]);
    end
    Ham = Hamiltonian(atoms, pspots, ecutwfc, options);
    return Ham
end
