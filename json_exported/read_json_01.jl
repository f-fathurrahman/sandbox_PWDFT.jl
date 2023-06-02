# Example reading JSON file

using JSON

# PAW radial integrator or PAWAtomicSphere

# For 1st species
json_data = JSON.parsefile("paw_radial_integrator_1.json")["paw_radial_integrator"]

json_nx = json_data["nx"]
json_ww = convert(Vector{Float64}, json_data["ww"])

shape_ylm = convert(Vector{Int64}, json_data["shape_ylm"])
json_ylm = reshape(convert(Vector{Float64}, json_data["ylm"]), shape_ylm...)

shape_wwylm = convert(Vector{Int64}, json_data["shape_wwylm"])
json_wwylm = reshape(convert(Vector{Float64}, json_data["wwylm"]), shape_wwylm...)


# pseudo_upf or PsPot_UPF
using JSON
json_data = JSON.parsefile("pseudo_upf_1.json")["pseudo_upf"];

shape_dion = convert(Vector{Int64}, json_data["shape_dion"]);
dion = reshape(convert(Vector{Float64}, json_data["dion"]), shape_dion...);

shape_qfuncl = convert(Vector{Int64}, json_data["shape_qfuncl"]);
qfuncl = reshape(convert(Vector{Float64}, json_data["qfuncl"]), shape_qfuncl...);

#
using JSON
json_data = JSON.parsefile("scf_mod.json")["scf_mod"];
rhoe_core = convert(Vector{Float64}, json_data["rho_core"]);

diff1 = Ham.rhoe_core  - rhoe_core;
maximum(abs.(diff1))

vltot = convert(Vector{Float64}, json_data["vltot"]);
