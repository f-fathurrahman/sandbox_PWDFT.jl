using PWDFT

include("read_pwscf_input_v2.jl")

println(ARGS[1])
atoms, meshk = read_pwscf_input(ARGS[1])

println(atoms)

println(meshk)