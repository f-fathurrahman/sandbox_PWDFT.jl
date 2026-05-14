mutable struct HamiltonianLAPW
    elk_input::ElkInput
    atoms::Atoms
    sym_vars::SymmetryVars
    specs_info::Vector{SpeciesInfo}
    atsp_vars::AtomicSpeciesVars
    mt_vars::MuffinTins
    apwlo_vars::APWLOVars
    sym_info::SymmetryInfo
    pw::PWGrid
    elec_chgst::ElectronicChargesStates
    gvec_full::GVectorsFull
    ffacg::Matrix{ComplexF64}
    cfunig::Vector{ComplexF64}
    cfunir::Vector{Float64}
    core_states::CoreStatesVars
    apwlo_ints::APWLOIntegrals
    nmat::Vector{Int64}
end

function HamiltonianLAPW()
    return
end
