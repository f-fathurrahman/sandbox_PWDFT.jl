using Infiltrator, Serialization

atoms = exfiltrated.atoms;
mt_vars = exfiltrated.mt_vars;
apwlo_vars = exfiltrated.apwlo_vars;
pw = exfiltrated.pw;
elec_chgst = exfiltrated.elec_chgst;

rhomt = exfiltrated.rhomt;
rhoir = exfiltrated.rhoir;
magmt = exfiltrated.magmt;
magir = exfiltrated.magir;
#serialize("rhomt_before.dat", rhomt)
#serialize("rhoir_before.dat", rhoir)
#serialize("magmt_before.dat", magmt)
#serialize("magir_before.dat", magir)

nstfv = elec_chgst.nstfv;
ispin = 1;
#ik = 1;

for ik in 1:pw.gvecw.kpoints.Nkpt

evecfv = deserialize("evecs_1st_ispin_$(ispin)_ik_$(ik).dat");
evecsv = deserialize("evecs_2nd_ispin_$(ispin)_ik_$(ik).dat");
apwalm = calc_match_coeffs(ik, atoms, pw, mt_vars, apwlo_vars);

#evecfv = deserialize("evecs_1st_ispin_$(ispin)_ik_$(ik)_elk.dat");
#apwalm_ia = deserialize("apwalm_ispin_$(ispin)_ik_$(ik)_elk.dat");

rhomagk!(
    ik, atoms, pw, mt_vars, apwlo_vars, elec_chgst,
    apwalm, evecfv, evecsv,
    rhomt, rhoir;
    magmt=magmt, magir=magir
)

end
