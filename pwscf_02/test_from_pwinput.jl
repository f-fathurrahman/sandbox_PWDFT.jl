using PWDFT

include("PWSCFInput.jl")

println(ARGS[1])
pwinput = PWSCFInput(ARGS[1])

Ham = Hamiltonian(
    pwinput.atoms, pwinput.pspfiles, pwinput.ecutwfc
)

psiks = rand_BlochWavefunc(Ham)
KS_solve_SCF!(Ham, psiks, mix_method="rpulay")
