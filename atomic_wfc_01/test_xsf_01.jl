filename = "PWINPUT_Al_v1";
#filename = "PWINPUT_H2O_oncv";
Ham, pwinput, Rhoe = PWDFT.prepare_Ham_from_pwinput(filename=filename);

# Prepare Haux
Nstates = Ham.electrons.Nstates;
Nspin = Ham.electrons.Nspin_wf;
Nkpt = Ham.pw.gvecw.kpoints.Nkpt;
Nkspin = Nkpt*Nspin;
Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin);
for ikspin in 1:Nkspin
    Haux[ikspin] = randn(ComplexF64, Nstates, Nstates);
    Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' );
end

# Prepare psiks
psiks = PWDFT.initwfc(Ham);
Nstates = Ham.electrons.Nstates;
ik = 1; # FIXED
psir = zeros(ComplexF64, prod(Ham.pw.Ns));
for ist in 1:Nstates
    psi = psiks[ik][:,ist];
    Ngwk = Ham.pw.gvecw.Ngw[ik];
    idx_gw2r = Ham.pw.gvecw.idx_gw2r;
    fill!(psir, 0.0)
    for igw in 1:Ngwk
        ip = idx_gw2r[ik][igw]
        psir[ip] = psi[igw]
    end
    G_to_R!(Ham.pw, psir)
    filexsf = "TEMP_wavefunc_ist_$ist.xsf"
    write_xsf(filexsf, Ham.atoms )
    write_xsf_data3d_crystal(filexsf, Ham.atoms, Ham.pw.Ns, real.(psir))
end
#electrons_Emin_Haux!(Ham, psiks=psiks, Haux=Haux, Rhoe=Rhoe)