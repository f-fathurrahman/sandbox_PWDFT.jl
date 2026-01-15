using LinearAlgebra
using PWDFT

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");
#Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv_nospinorb");

psiks = rand_BlochWavefunc(Ham);
Rhoe = Ham.rhoe; # alias
Rhoe[:,:], _ = atomic_rho_g( Ham,
    starting_magn = Ham.options.starting_magn,
    angle1 = Ham.options.angle1,
    angle2 = Ham.options.angle2 
);
update_from_rhoe!(Ham, psiks, Rhoe);

E_kin = calc_E_kin( Ham, psiks )
E_Ps_loc, E_Hartree, E_xc = calc_E_local( Ham, psiks )
E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )

# Compare with bra ket
K_psiks = op_K(Ham, psiks);
Vloc_psiks = op_V_loc(Ham, psiks);
Vpsnloc_psiks = op_V_Ps_nloc(Ham, psiks);

Nkpt = Ham.pw.gvecw.kpoints.Nkpt;
wk = Ham.pw.gvecw.kpoints.wk;
Focc = Ham.electrons.Focc;
Nstates = Ham.electrons.Nstates;

E_kin2 = 0.0;
for ik in 1:Nkpt, ist in 1:Nstates
    E_kin2 += wk[ik]*Focc[ist,ik]*real(dot(psiks[ik][:,ist], K_psiks[ik][:,ist]))
end

E_Ps_nloc2 = 0.0;
for ik in 1:Nkpt, ist in 1:Nstates
    E_Ps_nloc2 += wk[ik]*Focc[ist,ik]*real(dot(psiks[ik][:,ist], Vpsnloc_psiks[ik][:,ist]))
end