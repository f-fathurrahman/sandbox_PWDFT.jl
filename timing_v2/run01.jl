# To be copy-pasted

Ham_gth     = create_Ham_atom_Pt_gth();
Ham_oncv    = create_Ham_atom_Pt_oncv();
Ham_gbrv    = create_Ham_atom_Pt_gbrv();
Ham_paw_jth = create_Ham_atom_Pt_paw_jth();

Ham_gth     = create_Ham_O2_spinpol_gth();
Ham_oncv    = create_Ham_O2_spinpol_oncv();
Ham_gbrv    = create_Ham_O2_spinpol_gbrv();
Ham_paw_jth = create_Ham_O2_spinpol_paw_jth();


# PAWVariables
res = @be PAWVariables(Ham_paw_jth.atoms, Ham_paw_jth.pspots, Ham_paw_jth.pspotNL.nhm);

atoms = Ham_paw_jth.atoms;
pw = Ham_paw_jth.pw;
pspots = Ham_paw_jth.pspots;
#
pspotNL = Ham_paw_jth.pspotNL;
lmaxkb = pspotNL.lmaxkb;
nhm = pspotNL.nhm;
nh = pspotNL.nh;
indv = pspotNL.indv;
nhtolm = pspotNL.nhtolm
lpl = pspotNL.lpl;
lpx = pspotNL.lpx;
ap = pspotNL.ap;
res = @be PWDFT._prepare_aug_charges(
    atoms, pw, pspots, lmaxkb, nhm, nh, indv, nhtolm, lpl, lpx, ap
)