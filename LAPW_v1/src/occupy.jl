function occupy!(evalsv, occsv; NiterMax=1000)

    # XXX: automatic smearing width is disabled
    
    # find minimum and maximum eigenvalues
    e0 = evalsv[1,1]
    e1 = e0
    for ik in 1:nkpt
        for ist in 1:nstsv
            e = evalsv[ist,ik]
            if e < e0
                e0 = e
            end
            if e > e1
                e1 = e
            end
        end
    end
    
    if e0 < e0min
        println("WARNING: minimum eigenvalue less than minimum linearization energy : ", e0, e0min)
    end
    
    t1 = 1.0/swidth
    # determine the Fermi energy using the bisection method
    iterOcc = 0
    while iterOcc <= NiterMax
        efermi = 0.50*(e0 + e1)
        chg = 0.0
        for ik in 1:nkpt
            for ist in 1:nstsv
                e = evalsv[ist,ik]
                if e < e0min
                    occsv[ist,ik] = 0.0
                else
                    x = (efermi - e)*t1
                    occsv[ist,ik] = occmax*stheta_fd(x) # smearing
                    chg += wkpt[ik]*occsv[ist,ik]
                end 
            end
        end
        if chg < chgval
            e0 = efermi
        else
            e1 = efermi
        end
        if abs(e1-e0) < 1.e-12
            break
        end
        iterOcc += 1
    end
    if iterOcc > NiterMax 
        println("WARNING: could not find Fermi energy")
    end


    # find the density of states at the Fermi surface in units of
    # states/Hartree/unit cell
    fermidos = 0.0
    for ik in 1:nkpt
        for ist in 1:nstsv
            x = (evalsv[ist,ik] - efermi)*t1
            fermidos += wkpt[ik]*sdelta_fd(x)*t1
        end 
        if abs(occsv[nstsv,ik]) > epsocc
            println("WARNING: not enough empty states for k-point ", ik)
        end
    end
    fermidos = fermidos*occmax


#=
  ! estimate the indirect band gap (FC)
  e0=-1.d8
  e1=1.d8
  ikgap(1)=1
  ikgap(2)=1
  ! these loops are incorrectly ordered to fix a bug in versions 17 and 18 of the
  ! Intel compiler
  DO ist=1,nstsv
    DO ik=1,nkpt
      e=evalsv(ist,ik)
      IF(e.lt.efermi) THEN 
        IF(e.gt.e0) THEN 
          e0=e
          ikgap(1)=ik
        ENDIF 
      else
        IF(e.lt.e1) THEN 
          e1=e
          ikgap(2)=ik
        ENDIF 
      ENDIF 
    ENDDO 
  ENDDO 
  bandgap(1)=e1-e0

  ! write band gap to test file
  WRITE(*,*) 'estimated indirect band gap: ', bandgap(1)

  ! estimate the direct band gap
  e=1.d8
  ikgap(3)=1
  DO ik=1,nkpt
    e0=-1.d8
    e1=1.d8
    DO ist=1,nstsv
      t1=evalsv(ist,ik)
      IF(t1.le.efermi) THEN 
        IF(t1.gt.e0) e0=t1
      else
        IF(t1.lt.e1) e1=t1
      ENDIF 
    ENDDO 
    t1=e1-e0
    IF(t1.lt.e) THEN 
      e=t1
      ikgap(3)=ik
    ENDIF 
  ENDDO 
  bandgap(2)=e
=#

    return
end
