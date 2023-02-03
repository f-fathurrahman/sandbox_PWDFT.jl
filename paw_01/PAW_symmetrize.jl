
function PAW_symmetrize!( becsum )

    #REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin)
    #!! cross band occupations
    #!
    #! ... local variables
    #!
    #REAL(DP) :: becsym(nhm*(nhm+1)/2,nat,nspin) ! symmetrized becsum
    #REAL(DP) :: pref, usym, segno
    #REAL(DP) :: mb(3)
    #!
    #INTEGER :: ia,mykey,ia_s,ia_e 
    #!                       ! atoms counters and indexes
    #INTEGER :: is, nt       ! counters on spin, atom-type
    #INTEGER :: ma           ! atom symmetric to na
    #INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    #INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
    #           l_i, l_j, m_i, m_j
    #INTEGER :: m_o, m_u     ! counters for sums on m
    #INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    #INTEGER :: isym         ! counter for symmetry operation
    #INTEGER :: ipol, kpol

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
    
    D = Vector{Array{Float64,3}}(undef,4)
    D[1] = ones(1,1,48)
    D[2] = zeros(3,3,48)
    D[3] = zeros(5,5,48)
    D[4] = zeros(7,7,48)

    # IF( nsym==1 ) RETURN
    
    # d0(1,1,:) = 1._dp
    # D(0)%d => d0 ! d0(1,1,48)
    # D(1)%d => d1 ! d1(3,3,48)
    # D(2)%d => d2 ! d2(5,5,48)
    # D(3)%d => d3 ! d3(7,7,48)


# !
# ! => lm = l**2 + m
# ! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
# !       = lm + proj**2 + 2*l*proj
# !       = m + l**2 + proj**2 + 2*l*proj
# !        ^^^
# ! Known ih and m_i I can compute the index oh of a different m = m_o but
# ! the same augmentation channel (l_i = l_o, proj_i = proj_o):
# !  oh = ih - m_i + m_o
# ! this expression should be general inside pwscf.
# !
# !#define __DEBUG_PAW_SYM
# !
# !


    becsym = zeros(size(becsum))
    usym = 1.0/Nsyms

    for ispin in 1:Nspin, ia in 1:Natoms
        #
        isp = atm2species[ia]
        # No need to symmetrize non-PAW atoms
        # IF ( .NOT. upf(nt)%tpawp ) CYCLE
        if !pspot[isp].is_paw
            continue
        end
        #
        for ih in 1:nh[isp], jh in ih:nh[nt]
            # note: jh >= ih
            # ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ijh = ijtoh[ih,jh,isp]
            #
            lm_i = nhtolm[ih,isp]
            lm_j = nhtolm[jh,isp]
            #
            l_i = nhtol[ih,isp]
            l_j = nhtol[jh,isp]
            #
            m_i = lm_i - l_i^2
            m_j = lm_j - l_j^2
            #
            for isym in 1:Nsyms
                ma = irt[isym,ia]
                for m_o in 1:(2*l_i+1), m_u in (1:2*l_j+1)
                    oh = ih - m_i + m_o
                    uh = jh - m_j + m_u
                    ouh = ijtoh[oh,uh,nt]
                    # In becsum off-diagonal terms are multiplied by 2
                    # neutralize this factor and restore it later  
                    if oh == uh
                        pref = 2.0 * usym
                    else
                        pref = usym
                    end
                    #
                    becsym[ijh,ia,ispin += D[l_i+1][m_o,m_i,isym]*D[l_j+1][m_u,m_j,isym]*pref*becsum[ouh,ma,ispin]
                    # Starting index for D is 1
                end
            end # isym
            # Put the prefactor back in:  
            if ih == jh
                becsym[ijh,ia,ispin] *= 0.5
            end
        end # ih, jh
    end # ispin, ia

    # noncollinear spin is not yet implemented

    # Apply symmetrization:
    @views becsum[:,:,:] .= becsym[:,:,:]

    return

end