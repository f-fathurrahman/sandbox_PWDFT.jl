using Revise, PWDFT
using LinearAlgebra

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..");
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials");

includet("create_Ham_various.jl")
