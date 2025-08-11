include("setup_path.jl")
include("test_cg_01.jl")

# choose a system
Ham = create_Ham_Fe_bcc_smearing_paw_jth();

# Our current direct minimization algorithm need quite good
# starting point psiks.
# We start with psiks from SCF
psiks = rand_BlochWavefunc(Ham);
PWDFT.electrons_scf_G!(Ham, psiks=psiks, NiterMax=2, betamix=0.1)
# psiks will be updated after this call
# set NiterMax to 1 or 2

# for some systems electrons_scf_G! is more stable


# Now run the direct minimization algorithm
main_cg_01(Ham, psiks=psiks)
# sometimes it stuck to local min, so it is useful to run it again.

# compare with SCF
PWDFT.electrons_scf_G!(Ham, psiks=psiks)

