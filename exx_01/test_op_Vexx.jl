using Printf
using LinearAlgebra: norm, inv, dot
using FFTW
using Serialization: serialize, deserialize

includet("cryst_to_cart.jl")
includet("exx_grid_check.jl")
includet("exx_qgrid_init.jl")
includet("scale_sym_ops.jl")
includet("rotate_grid_point.jl")
includet("exx_set_symm.jl")
includet("calc_exx_divergence.jl")
includet("EXXVariables.jl")
includet("g2_convolution.jl")
includet("op_Vexx.jl")

function debug_op_Vexx()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    exx = EXXVariables(Ham, pwinput)
    psiks = deserialize("psiks_nox_noc.jldat")
    set_exx_buffer!(Ham, exx, psiks)

    psiks_in = psiks # rand_BlochWavefunc(Ham)
    psiks_out = zeros_BlochWavefunc(Ham)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ham.ispin = 1
    for ik in 1:Nkpt
        Ham.ik = ik
        psi = psiks_in[ik]
        Vpsi = psiks_out[ik]
        op_Vexx!(Ham, exx, psi, Vpsi)
    end

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    wk = Ham.pw.gvecw.kpoints.wk
    ene_exx = 0.0
    for ik in 1:Nkpt
        psi = psiks_in[ik]
        Vpsi = psiks_out[ik]
        for ist in 1:Nstates
            ene_exx += wk[ik] * Focc[ist,ik] * dot(psi[:,ist], Vpsi[:,ist])
        end
    end
    println("ene_exx (in Ry) = ", 2*ene_exx)

    @infiltrate
end
