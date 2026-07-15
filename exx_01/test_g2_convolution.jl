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
includet("calc_exx_divergence.jl")
includet("EXXVariables.jl")
includet("g2_convolution.jl")

function debug_g2_convolution()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    exx = EXXVariables(Ham, pwinput)
    psiks = deserialize("psiks_nox_noc.jldat")
    set_exx_buffer!(Ham, exx, psiks)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    nqs = exx.nqs
    Ng = exx.gvec.Ng
    fac = zeros(Float64, Ng)

    g2_convolution!(exx, Ham.pw.LatVecs, Ham.pw.gvecw.kpoints.k[:,1], exx.xkq[:,1], fac)
    println("sum fac (in Ry) = ", sum(fac)*2.0)

    #for ik in 1:Nkpt
    #    fill!(fac, 0.0)
    #    g2_convolution!(exx, Ham.pw.LatVecs, Ham.pw.gvecw.kpoints.k[:,ik], exx.xkq[:,1], fac)
    #    println("sum(fac) = ", sum(fac))
    #end

    @infiltrate
end
