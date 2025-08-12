# To be copy-pasted in REPL

using Revise, Infiltrator
using LinearAlgebra
using Printf
using Serialization
using Random
using PWDFT

includet("prepare_Ham_various.jl")

includet("smearing.jl")
includet("occupations.jl")

includet("electrons_Emin_Haux.jl")

# choose a system
#Ham = create_Ham_Fe_bcc_smearing_paw_jth();
Ham = create_Ham_Pt_fcc_smearing();

# Our current direct minimization algorithm need quite good
# starting point psiks.
# We start with psiks from SCF
psiks = rand_BlochWavefunc(Ham);
PWDFT.electrons_scf_G!(Ham, psiks=psiks, NiterMax=2, betamix=0.1)
# psiks will be updated after this call
# set NiterMax to 1 or 2

# for some systems electrons_scf_G! is more stable

# Now run the direct minimization algorithm
electrons_Emin_Haux!(Ham, psiks=psiks)
# sometimes it stuck to local min, so it is useful to run it again.

# compare with SCF
PWDFT.electrons_scf_G!(Ham, psiks=psiks)

