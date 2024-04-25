# will modify some fields of atoms
function findsymcrys!(
    atoms;
    tshift=true,
    epslat=1e-6,
    symtype=1
)

    #=
    IMPLICIT NONE 
    ! local variables
    INTEGER :: ia,ja,is,js
    INTEGER :: isym,nsym,i,n
    INTEGER :: lspl(48),lspn(48),ilspl
    REAL(8) :: v0(3),v1(3),v2(3),t1
    REAL(8) :: apl(3,maxatoms,maxspecies)
    ! ALLOCATABLE arrays
    INTEGER, ALLOCATABLE :: iea(:,:,:)
    REAL(8), ALLOCATABLE :: vtl(:,:)
    =#

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    LatVecs = atoms.LatVecs
    atm2species = atoms.atm2species
    atposc = atoms.positions
    atposl = inv(atoms.LatVecs)*atposc

    MAX_SYM_CRYS = 172

    ieqatom = zeros(Int64, Natoms, MAX_SYM_CRYS)    
    eqatoms = zeros(Bool, Natoms, Natoms)

    # store position of first atom
    v0 = zeros(Float64, 3)
    if Natoms > 0
        @views v0[:] = atposl[:,1]
    end

    # find the smallest set of atoms
    # XXX: Need this?
    Natoms_per_species = zeros(Int64, Nspecies)
    for isp in 1:Nspecies
        Natoms_per_species[isp] = count(atm2species .== isp)
    end
    isp_smallest = 1
    for isp in 1:Nspecies
        if Natoms_per_species[isp] < Natoms_per_species[isp_smallest]
            isp_smallest = isp
        end
    end

    v1 = zeros(Float64, 3)
    # atomis.positions will be shifted
    if tshift
        # shift basis so that the first atom in the smallest atom set is at the origin
        #
        # find the first atom with the species of smallest number of atoms
        for ia in Natoms
            isp = atm2species[ia]
            if isp == isp_smallest
                v1[:] = atposl[:,ia]
                break
            end
        end
        #
        for ia in 1:Natoms
            # shift atom
            atposl[:,ia] = atposl[:,ia] - v1[:]
            # map lattice coordinates back to [0,1)
            r3frac!(atposl[:,ia])
            # determine the new Cartesian coordinates
            atposc[:,ia] = LatVecs * atposl[:,ia]
        end 
    end

    # determine possible translation vectors from smallest set of atoms

    #n = max( natoms(is)*natoms(is), 1 )
    vtl = zeros(Float64, 3, Natoms_per_species[isp_smallest]*Natoms_per_species[isp_smallest])
    n = 1
    vtl[:,1] .= 0.0 # not really needed
    for ia in 1:Natoms_per_species[isp_smallest]
        for ja in 2:Natoms_per_species[isp_smallest]
            # compute difference between two atom vectors
            v1[:] = atposl[:,ia] - atposl[:,ja]
            # map lattice coordinates to [0,1)
            r3frac!(v1)
            # 
            # .... tefield specific treatment is removed
            # 
            for i in 1:n
                t1 = abs( vtl[1,i] - v1[1] ) + abs( vtl[2,i] - v1[2] ) + abs( vtl[3,i] - v1[3] )
                if t1 < epslat
                    @goto LABEL10
                end
            end 
            n += 1
            vtl[:,n] = v1[:]
            @label LABEL10
        end
    end


    # no translations required when symtype=0,2 (F. Cricchio)
    if symtype != 1
        n = 1
    end

    println("n = ", n)
    eqatoms = zeros(Bool, Natoms, Natoms)
    apl = zeros(Float64, 3, Natoms)
    nsymcrys = 0
    lspl = zeros(Int64, 48)
    lspn = zeros(Int64, 48)
    iea = zeros(Int64, Natoms, 48)

    # loop over all possible translations
    for i in 1:n
        # construct new array with translated positions
        for ia in 1:Nspecies
            apl[:,ia] = atposl + vtl[:,i]
        end
        # find the symmetries for current translation
        nsym = findsym!(atposl, apl, nsym, lspl, lspn, iea)
        for isym in 1:nsym
            nsymcrys = nsymcrys + 1
            if nsymcrys > MAX_SYM_CRYS 
                error("Too many nsymcrys = $(nsymcrys)")
            end
            @views vtlsymc[:,nsymcrys] = vtl[:,i]
            lsplsymc[nsymcrys] = lspl[isym]
            lspnsymc[nsymcrys] = lspn[isym]
            for ia in 1:Natoms
                ja = iea[ia,isym]
                ieqatom[ia,nsymcrys] = ja
                eqatoms[ia,ja] = true
                eqatoms[ja,ia] = true
            end 
        end
    end

    tsyminv = false
    for isym in 1:nsymcrys
        i = lsplsymc[isym]
        # check if inversion symmetry is present
        if symlat[i] == -symlat[1]
            tsyminv = true
            # make inversion the second symmetry element (the identity is the first)
            v1[:] = vtlsymc[:,isym]; vtlsymc[:,isym] = vtlsymc[:,2]; vtlsymc[:,2] = v1(:)
            i = lsplsymc[isym]; lsplsymc[isym] = lsplsymc[2]; lsplsymc[2] = i
            i = lspnsymc[isym]; lspnsymc[isym] = lspnsymc[2]; lspnsymc[2] = i
            for ia in 1:Natoms
                i = ieqatom[ia,isym]
                ieqatom[ia,isym] = ieqatom[ia,2]
                ieqatom[ia,2] = i
            end
            break
        end
    end


    # if inversion exists THEN  shift basis so that inversion center is at origin
    if tsyminv && tshift 
        v1[:] = v1[:]/2.0
        for ia in 1:Natoms
            # shift atom
            atposl(:,ia,is) = atposl(:,ia,is) + v1(:)
        ! map lattice coordinates back to [0,1)
        CALL r3frac(epslat,atposl(:,ia,is))
        ! map lattice coordinates to [-0.5,0.5)
        DO i = 1,3
          IF( atposl(i,ia,is) > 0.5d0 ) atposl(i,ia,is) = atposl(i,ia,is) - 1.d0
        ENDDO 
        ! determine the new Cartesian coordinates
        CALL r3mv(avec, atposl(:,ia,is), atposc(:,ia,is))
      ENDDO 
    ENDDO 
    ! recalculate crystal symmetry translation vectors
    DO isym = 1,nsymcrys
      ilspl = isymlat(lsplsymc(isym))
      v2(:) = symlat(:,1,ilspl)*v1(1) + symlat(:,2,ilspl)*v1(2) + symlat(:,3,ilspl)*v1(3)
      vtlsymc(:,isym) = vtlsymc(:,isym) - v1(:) + v2(:)
      CALL r3frac(epslat, vtlsymc(:,isym))
    ENDDO 
  ENDIF 

  ! translation vector in Cartesian coordinates
  DO isym = 1,nsymcrys
    CALL r3mv(avec, vtlsymc(:,isym), vtcsymc(:,isym))
  ENDDO 

  ! set flag for zero translation vector
  DO isym = 1,nsymcrys
    t1 = abs(vtlsymc(1,isym)) + abs(vtlsymc(2,isym)) + abs(vtlsymc(3,isym))
    IF(t1 < epslat) THEN 
      tv0symc(isym) = .true.
    else
      tv0symc(isym) = .false.
    ENDIF 
  ENDDO 

  ! check inversion does not include a translation
  IF(tsyminv) THEN 
    IF(.not. tv0symc(2) ) tsyminv = .false.
  ENDIF 

  IF( natmtot > 0 ) THEN 
    v1(:) = atposl(:,1,1) - v0(:)
    t1 = abs(v1(1)) + abs(v1(2)) + abs(v1(3))
    if( t1 > epslat ) THEN 
      WRITE(*,*)
      WRITE(*,'("Info(findsymcrys): atomic basis shift (lattice) :")')
      WRITE(*,'(3G18.10)') v1(:)
      WRITE(*,'("See GEOMETRY.OUT for new atomic positions")')
    ENDIF 
  ENDIF 

    return

#=
=#


end
