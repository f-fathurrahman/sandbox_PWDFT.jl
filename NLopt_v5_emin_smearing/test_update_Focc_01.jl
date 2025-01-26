using LinearAlgebra
using Printf

using PWDFT

include("smearing.jl")
include("occupations.jl")

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");

# This will take into account whether the overlap operator is needed or not
psiks = rand_BlochWavefunc(Ham)

use_smearing = false
kT = 0.0
if pwinput.occupations == "smearing"
    use_smearing = true
    kT = pwinput.degauss*0.5 # convert from Ry to Ha
end

if pwinput.nspin == 2
    starting_magnetization = pwinput.starting_magnetization
else
    starting_magnetization = nothing
end

Nspin = Ham.electrons.Nspin
Nkpt = Ham.pw.gvecw.kpoints.Nkpt
Nkspin = Nkpt*Nspin
Nstates = Ham.electrons.Nstates
ebands = Ham.electrons.ebands
Focc = Ham.electrons.Focc
wk = Ham.pw.gvecw.kpoints.wk
Nelectrons = Ham.electrons.Nelectrons


Haux = Vector{Matrix{Float64}}(undef, Nkspin)
for ikspin in 1:Nkspin
    Haux[ikspin] = randn(Float64, Nstates, Nstates)
    Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
    ebands[:,ikspin] = eigvals(Haux[ikspin])
end

update_Focc!(
    Focc,
    smear_fermi, smear_fermi_entropy,
    ebands,
    Nelectrons,
    kT,
    Nkpt,
    wk
)