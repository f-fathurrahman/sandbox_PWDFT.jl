using PWDFT

function main()
    LatVecs = [16.0 0.0 0.0
               0.0 16.0 0.0
               0.0 0.0 16.0] 

    gmaxvr = 12.0
    ecutrho = 0.5*gmaxvr^2
    
    gkmax = 4.2452830189
    ecutwfc = 0.5*gkmax^2

    dual = ecutrho/ecutwfc
    println("dual = ", dual)

    pw = PWGrid(ecutwfc, LatVecs, dual=dual)
    println(pw)

    #ecutwfc = 0.5*gkmax^2
    #kpoints = KPoints( 1, (1,1,1), zeros(3,1), [1.0], pw.RecVecs )
    #gvecw = PWDFT.init_gvecw( ecutwfc, pw.gvec, kpoints )
    #println(gvecw.Ngw)
end

main()