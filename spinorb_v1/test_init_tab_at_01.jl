using Revise
includet("init_tab_at.jl")

function test_main(; filename=nothing, do_export_data=false)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)
    tab_at = init_tab_at(Ham.pspots[1], Ham.pw)
end
#=
Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT_oncv");
=#