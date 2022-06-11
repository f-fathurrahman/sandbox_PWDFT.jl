using PWDFT

include("PWSCFInput.jl")

println(ARGS[1])
pwinput = PWSCFInput(ARGS[1])

println(pwinput)