using OffsetArrays

function test1()
    lmaxo = 2
    Nspecies = 2
    nrcmtmax = 10
    rlmt = OffsetArray( zeros(nrcmtmax, 2*(lmaxo+2), Nspecies), 
        1:nrcmtmax, -lmaxo-1:lmaxo+2, 1:Nspecies
    )
    rlmt2 = 2*rlmt
    println(typeof(rlmt))
    #println(typeof(rlmt2))
end
@time test1()
#@time test1()

function test2()
    lmaxo = 2
    Nspecies = 2
    nrcmt = zeros(Int64,Nspecies)
    nrcmt[1] = 10
    nrcmt[2] = 11
    rlmt = Array{OffsetMatrix{Float64, Matrix{Float64}},1}(undef,Nspecies)

    for isp in 1:Nspecies
        rlmt[isp] = OffsetArray( zeros(nrcmt[isp], 2*(lmaxo+2)), 
            1:nrcmt[isp], -lmaxo-1:lmaxo+2
        )
    end
    println(typeof(rlmt[1]))
end
@time test2()