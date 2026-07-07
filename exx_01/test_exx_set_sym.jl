using Infiltrator

includet("scale_sym_ops.jl")
includet("rotate_grid_point.jl")
includet("exx_set_symm.jl")

function debug_main()
    filename = "PWINPUT_AlAs"
    Ham, pwinput = init_Ham_from_pwinput( filename = filename )

    rir = exx_set_symm( Ham.sym_info, Ham.pw.Ns )

    @infiltrate
end
