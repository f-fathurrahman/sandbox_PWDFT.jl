using LAPWDFT

function main()
    lmaxi = 1
    lmaxo = 6
    shtmat = SphericalHarmonicTransform(lmaxi, lmaxo)
    display(shtmat.rbshto); println()
    display(shtmat.rbshti); println()
end

main()