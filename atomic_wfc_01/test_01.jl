using Revise, PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..");
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials");

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

#=
Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");
=#

includet("init_tab_at.jl")
includet("atomic_wfc_01.jl")

function test_main(; filename=nothing, do_export_data=false)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    tab_at = init_tab_at(Ham.pspots[1], Ham.pw)
end


Nstates = Ham.electrons.Nstates;
Nspin = Ham.electrons.Nspin_wf;
Nkpt = Ham.pw.gvecw.kpoints.Nkpt;
psiks = zeros_BlochWavefunc(Ham);
Natomwfc = calc_Natomwfc(Ham.atoms, Ham.pspots);
println("Natomwfc = ", Natomwfc);
println("Nstates = ", Nstates);
for ispin in 1:Nspin, ik in 1:Nkpt
    ikspin = ik + (ispin-1)*Nkpt
    @views atomic_wfc!(ik, Ham.atoms, Ham.pspots, Ham.pw, psiks[ikspin][:,1:Natomwfc]);
end

for ispin in 1:Nspin, ik in 1:Nkpt
    ikspin = ik + (ispin-1)*Nkpt
    ortho_sqrt!(Ham, psiks[ikspin])
end

Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin);
for ikspin in 1:Nkspin
    Haux[ikspin] = randn(ComplexF64, Nstates, Nstates);
    Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' );
end


"""

Ham.rhoe[:,:], _ = atomic_rho_g(
    Ham,
    starting_magn = Ham.options.starting_magn,
    angle1 = Ham.options.angle1,
    angle2 = Ham.options.angle2
)
update_from_rhoe!(Ham, psiks, Ham.rhoe)

ik = 1; ispin = 1;
for ispin in 1:Nspin, ik in 1:Nkpt
    Ham.ispin = ispin
    Ham.ik = ik
    ikspin = ik + (ispin-1)*Nkpt
    psi = psiks[ikspin];
    Hpsi = op_H(Ham, psiks[ikspin]);
    Hsub = psi' * Hpsi;
    λ, U = eigen(Hsub);
end
"""