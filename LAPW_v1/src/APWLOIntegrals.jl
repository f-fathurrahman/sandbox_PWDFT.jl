mutable struct APWLOIntegrals
    haa::Vector{OffsetArray{Float64, 5, Array{Float64, 5}}}
    hloa::Vector{OffsetArray{Float64, 4, Array{Float64, 4}}}
    hlolo::Vector{Array{Float64,3}}
    oalo::Vector{Matrix{Float64}}
    ololo::Vector{Matrix{Float64}}
end

function APWLOIntegrals( atoms, mt_vars, apwlo_vars )
    Natoms = atoms.Natoms
    nlomax = apwlo_vars.nlomax
    apwordmax = apwlo_vars.apwordmax
    lmaxapw = mt_vars.lmaxapw
    lmmaxo = mt_vars.lmmaxo

    oalo = Vector{Matrix{Float64}}(undef, Natoms)
    for ia in 1:Natoms
        oalo[ia] = zeros(Float64, apwordmax, nlomax)
    end

    ololo = Vector{Matrix{Float64}}(undef, Natoms)
    for ia in 1:Natoms
        ololo[ia] = zeros(Float64, nlomax, nlomax)
    end

    haa = Vector{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, Natoms)
    for ia in 1:Natoms
        haa[ia] = OffsetArray(
            zeros(Float64, lmmaxo, apwordmax, lmaxapw+1, apwordmax, lmaxapw+1),
            1:lmmaxo, 1:apwordmax, 0:lmaxapw, 1:apwordmax, 0:lmaxapw
        )
    end

    hloa = Vector{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, Natoms)
    for ia in 1:Natoms
        hloa[ia] = OffsetArray(
            zeros(Float64, lmmaxo, apwordmax, lmaxapw+1, nlomax),
            1:lmmaxo, 1:apwordmax, 0:lmaxapw, 1:nlomax
        )
    end

    hlolo = Vector{Array{Float64,3}}(undef, Natoms)
    for ia in 1:Natoms
        hlolo[ia] = zeros(Float64, lmmaxo, nlomax, nlomax)
    end

    return APWLOIntegrals( haa, hloa, hlolo, oalo, ololo )
end

function calc_apwlo_integrals!(atoms, mt_vars, apwlo_vars, vsmt, apwlo_ints)
    hmlrad!( atoms, mt_vars, apwlo_vars, vsmt, apwlo_ints.haa, apwlo_ints.hloa, apwlo_ints.hlolo )
    olprad!( atoms, mt_vars, apwlo_vars, apwlo_ints.oalo, apwlo_ints.ololo )
    return
end
