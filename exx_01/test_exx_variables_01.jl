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

function debug_main()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    exx = EXXVariables(Ham, pwinput)
    psiks = deserialize("psiks_nox_noc.jldat")
    set_exx_buffer!(Ham, exx, psiks)

    # Display q-points
    q = zeros(Float64, 3)
    xkq = zeros(Float64, 3)
    xk = Ham.pw.gvecw.kpoints.k
    nqs = exx.nqs
    current_ik = 1
    for iq in 1:nqs
        ikq = exx.index_xkq[current_ik,iq]
        ik = exx.index_xk[ikq]
        xkq[:] = exx.xkq[:,ikq]
        q[:] = xk[:,current_ik] - xkq[:]
        println("q = ", q)
    end

    @infiltrate
end
