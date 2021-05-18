function genylmg!(lmaxo, G, ylmg)
    Ng = size(G,2)
    for ig in 1:Ng
        @views genylmv!(lmaxo, G[:,ig], ylmg[:,ig] )
    end
    return
end