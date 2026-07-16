using Printf
using LinearAlgebra: norm, inv, dot
using FFTW
using Serialization: serialize, deserialize

function debug_diag_exx()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    
    Ham.exx = EXXVariables(Ham, pwinput)
    psiks = deserialize("psiks_nox_noc.jldat")
    #psiks = deserialize("psiks_nox_noc_v01.jldat")
    Ham.exx.is_active = true
    set_exx_buffer!(Ham, psiks)

    psiks_in = psiks # rand_BlochWavefunc(Ham)
    psiks_out = zeros_BlochWavefunc(Ham)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ham.ispin = 1
    for ik in 1:Nkpt
        Ham.ik = ik
        psi = psiks_in[ik]
        Hpsi = psiks_out[ik]
        op_H!(Ham, psi, Hpsi)
    end

    # Test calculating energy
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    wk = Ham.pw.gvecw.kpoints.wk
    ene = 0.0
    for ik in 1:Nkpt
        psi = psiks_in[ik]
        Hpsi = psiks_out[ik]
        for ist in 1:Nstates
            ene += wk[ik] * Focc[ist,ik] * dot(psi[:,ist], Hpsi[:,ist])
        end
    end
    println("ene (in Ry) = ", 2*ene)

    @infiltrate
end
