function linmin_grad!(
    Ham::Hamiltonian, psiks::BlochWavefunc,
    g::BlochWavefunc, d::BlochWavefunc,
    psic::BlochWavefunc, gt::BlochWavefunc;
    αt = 3e-5
)
    Nkspin = length(psiks)
    for i in 1:Nkspin
        psic[i] = psiks[i] + αt*d[i]
        ortho_sqrt!( psic[i] )
    end

    Rhoe = calc_rhoe(Ham, psic)
    update!(Ham, Rhoe)

    calc_grad!( Ham, psic, gt )

    denum = 2.0*real( dot(g .- gt, d) )
    if denum != 0.0
        α = abs( αt * 2.0*real( dot(g, d) )/denum )
    else
        α = 0.0
    end

    return α

end