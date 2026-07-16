using Printf
using LinearAlgebra: norm, inv
using FFTW
using Serialization: serialize, deserialize

includet("cryst_to_cart.jl")
includet("exx_grid_check.jl")
includet("exx_qgrid_init.jl")
includet("scale_sym_ops.jl")
includet("rotate_grid_point.jl")
includet("exx_set_symm.jl")
includet("EXXVariables.jl")
includet("calc_exx_divergence.jl")

function debug_exx_divergence()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    exx = EXXVariables(Ham, pwinput)
    psiks = deserialize("psiks_nox_noc.jldat")
    set_exx_buffer!(Ham, exx, psiks)

    exxdiv = calc_exx_divergence(
        pw, exx.Nq1, exx.Nq2, exx.Nq3,
        exx.x_gamma_extrapolation,
        exx.erf_scrlen, exx.erfc_scrlen, exx.yukawa, exx.grid_factor, exx.nqs,
        use_regularization = exx.use_regularization
    )

    @infiltrate
end
