function rf_unpack!()
    # unpack the function
    for ia in 1:Natoms
        isp = atm2species[ia]
        #call dcopy(np(is),v(n+1),1,rfmt(:,ias),1)
        n += np[isp]
    end
    #call dcopy(ngtot,v(n+1),1,rfir,1)
    n += ngtot
    return n
end
