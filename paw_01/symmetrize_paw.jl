function symmetrize_paw( becsum )

    # USE lsda_mod,          ONLY : nspin
    # USE cell_base,         ONLY : at, bg
    # USE noncollin_module,  ONLY : nspin_lsda
    # USE spin_orb,          ONLY : domag
    # USE uspp_param,        ONLY : nhm
    # USE ions_base,         ONLY : nat, ityp
    # USE symm_base,         ONLY : nsym, irt, d1, d2, d3, t_rev, sname, s, &
    #                               invs, inverse_s
    # USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    # USE uspp_param,        ONLY : nh, upf
    # USE io_global,         ONLY : stdout, ionode
    
    # REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin)
    
    # !! cross band occupations
    # !
    # ! ... local variables
    # !
    # REAL(DP) :: becsym(nhm*(nhm+1)/2,nat,nspin) ! symmetrized becsum
    # REAL(DP) :: pref, usym, segno
    # REAL(DP) :: mb(3)
    
    !
    INTEGER :: ia,mykey,ia_s,ia_e 
    !                       ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym         ! counter for symmetry operation
    INTEGER :: ipol, kpol
    INTEGER :: table(48,48)

    # !
    # ! The following mess is necessary because the symmetrization operation
    # ! in LDA+U code is simpler than in PAW, so the required quantities are
    # ! represented in a simple but not general way.
    # ! I will fix this when everything works.
    # REAL(DP), TARGET :: d0(1,1,48)
    # TYPE symmetrization_tensor
    #     REAL(DP),POINTER :: d(:,:,:)
    # END TYPE symmetrization_tensor
    # TYPE(symmetrization_tensor) :: D(0:3)
    
    IF( nsym==1 ) RETURN

    d0(1,1,:) = 1.0
    D(0)%d => d0 # d0(1,1,48)
    D(1)%d => d1 # d1(3,3,48)
    D(2)%d => d2 # d2(5,5,48)
    D(3)%d => d3 # d3(7,7,48)

    #!
    #! => lm = l**2 + m
    #! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
    #!       = lm + proj**2 + 2*l*proj
    #!       = m + l**2 + proj**2 + 2*l*proj
    #!        ^^^
    #! Known ih and m_i I can compute the index oh of a different m = m_o but
    #! the same augmentation channel (l_i = l_o, proj_i = proj_o):
    #!  oh = ih - m_i + m_o
    #! this expression should be general inside pwscf.
    #!
    #!#define __DEBUG_PAW_SYM
    #!
    #!

    becsym(:,:,:) = 0._dp
    usym = 1._dp / DBLE(nsym)

    !
    DO is = 1, nspin_lsda
       !
       atoms: DO ia = ia_s, ia_e
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .NOT. upf(nt)%tpawp ) CYCLE
          !
          DO ih = 1, nh(nt)
          DO jh = ih, nh(nt) ! note: jh >= ih
             !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
             ijh = ijtoh(ih,jh,nt)
             !
             lm_i  = nhtolm(ih,nt)  
             lm_j  = nhtolm(jh,nt)  
             !  
             l_i   = nhtol(ih,nt)  
             l_j   = nhtol(jh,nt)  
             !  
             m_i   = lm_i - l_i**2  
             m_j   = lm_j - l_j**2  
             !  
             DO isym = 1,nsym  
                ma = irt(isym,ia)  
                DO m_o = 1, 2*l_i +1  
                DO m_u = 1, 2*l_j +1  
                   oh = ih - m_i + m_o  
                   uh = jh - m_j + m_u  
                   ouh = ijtoh(oh,uh,nt)  
                   ! In becsum off-diagonal terms are multiplied by 2, I have  
                   ! to neutralize this factor and restore it later  
                   IF ( oh == uh ) THEN  
                      pref = 2._dp * usym  
                   ELSE  
                      pref = usym  
                   ENDIF  
                   !  
                   becsym(ijh, ia, is) = becsym(ijh, ia, is) &  
                       + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &  
                         * pref * becsum(ouh, ma, is)  
                ENDDO ! m_o  
                ENDDO ! m_u  
             ENDDO ! isym  
             ! Put the prefactor back in:  
             IF ( ih == jh ) becsym(ijh,ia,is) = .5_dp * becsym(ijh,ia,is)
             !
          ENDDO ! ih
          ENDDO ! jh
          !
       ENDDO atoms ! nat
       !
    ENDDO ! nspin
    !
    !
    IF (nspin==4 .AND. domag) THEN
       !
       CALL inverse_s( )
       !
       becsym(:,:,2:4) = 0._dp
       DO ia = 1, nat
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .NOT. upf(nt)%tpawp ) CYCLE
          !
          !  Bring the magnetization in the basis of the crystal
          !        
          DO ijh=1,(nh(nt)*(nh(nt)+1))/2
             DO ipol=1,3
                mb(ipol) = becsum(ijh,ia,ipol+1)
             ENDDO
             DO ipol=1,3
                becsum(ijh,ia,ipol+1) = bg(1,ipol)*mb(1)+bg(2,ipol)*mb(2) + &
                                        bg(3,ipol)*mb(3) 
             ENDDO
          ENDDO
       ENDDO
       !
       atoms_1: DO ia = ia_s, ia_e
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .NOT. upf(nt)%tpawp ) CYCLE
          !
          DO ih = 1, nh(nt)
          DO jh = ih, nh(nt) ! note: jh >= ih
             !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
             ijh = ijtoh(ih,jh,nt)
             !
             lm_i  = nhtolm(ih,nt)
             lm_j  = nhtolm(jh,nt)  
             !  
             l_i   = nhtol(ih,nt)  
             l_j   = nhtol(jh,nt)  
             !  
             m_i   = lm_i - l_i**2  
             m_j   = lm_j - l_j**2  
             !  
             DO isym = 1, nsym  
                ma = irt(isym,ia)  
                DO m_o = 1, 2*l_i +1  
                DO m_u = 1, 2*l_j +1  
                   oh = ih - m_i + m_o  
                   uh = jh - m_j + m_u  
                   ouh = ijtoh(oh,uh,nt)
                   ! In becsum off-diagonal terms are multiplied by 2, I have  
                   ! to neutralize this factor and restore it later  
                   IF ( oh == uh ) THEN  
                       pref = 2._dp * usym  
                   ELSE  
                       pref = usym  
                   ENDIF  
                   !  
                   segno=1.0_DP  
                   IF (sname(isym)(1:3)=='inv') segno = -segno  
                   IF (t_rev(isym)==1)  segno = -segno  
                   !  
                   DO is = 1, 3  
                   DO kpol = 1, 3  
                      becsym(ijh,ia,is+1) = becsym(ijh,ia,is+1)            &  
                         + D(l_i)%d(m_o,m_i,isym) * D(l_j)%d(m_u,m_j,isym) &  
                           * pref * becsum(ouh,ma,kpol+1)                  &  
                           * s(kpol,is,invs(isym)) * segno  
                   ENDDO  
                   ENDDO  
                ENDDO ! m_o  
                ENDDO ! m_u  
             ENDDO ! isym  
             !  
             ! Put the prefactor back in:  
             IF ( ih == jh ) becsym(ijh,ia,2:4) = .5_dp * becsym(ijh,ia,2:4)
             !
          ENDDO ! ih
          ENDDO ! jh
       ENDDO atoms_1 ! nat
       !
       DO ia = ia_s, ia_e
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .NOT. upf(nt)%tpawp ) CYCLE
          !
          !  ... Bring the magnetization in cartesian basis
          !        
          DO ijh = 1, (nh(nt)*(nh(nt)+1))/2
             DO ipol = 1, 3
                mb(ipol) = becsym(ijh,ia,ipol+1)
             ENDDO
             DO ipol = 1, 3
                becsym(ijh,ia,ipol+1) = at(ipol,1)*mb(1)+at(ipol,2)*mb(2) + &
                                        at(ipol,3)*mb(3)
             ENDDO
          ENDDO
       ENDDO
       !
    ENDIF
    !
#if defined(__MPI)
    IF ( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum( becsym, intra_image_comm )
#endif
    !
#if defined(__DEBUG_PAW_SYM)
    WRITE (stdout,*) "------------"
    IF (ionode) THEN
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin
            WRITE (*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            WRITE (stdout,"(1f10.3)", ADVANCE='no') becsym(ijh,ia,is)
        ENDDO
            WRITE (stdout,*)
        ENDDO
            WRITE (stdout,*)
        ENDDO
    ENDIF
    WRITE (stdout,*) "------------"
#endif
    !
    ! Apply symmetrization:
    becsum(:,:,:) = becsym(:,:,:)
    !
    CALL stop_clock( 'PAW_symme' )
    !
END SUBROUTINE PAW_symmetrize