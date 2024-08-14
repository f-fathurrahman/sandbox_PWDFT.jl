#=
!INPUT/OUTPUT PARAMETERS:
!   thr    : .true. if the matrix h is real valued (in,logical)
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients (in,complex(ngkmax,apwordmax,lmmaxapw))
!   ld     : leading dimension of h (in,integer)
!   h      : Hamiltonian matrix (inout,complex(*))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
=#
function hmlaa!(ia, atoms, pw, mt_vars, apwlo_vars, apwalm, h)
    #=
    USE m_atoms, ONLY: idxis
  USE m_gkvectors, ONLY: ngkmax
  USE m_apwlo, ONLY: apwordmax, apwfr, apwdfr, apword, lmoapw
  USE m_muffin_tins, ONLY: nrmt, idxlm, rmt, lmaxo, lmaxapw, lmmaxapw
  USE m_hamiltonian, ONLY: haa, gntyry
  IMPLICIT NONE 
  ! arguments
  LOGICAL, INTENT(in) :: thr
  INTEGER, INTENT(in) :: ias,ngp
  COMPLEX(8), INTENT(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
  INTEGER, INTENT(in) :: ld
  COMPLEX(8), INTENT(inout) :: h(*)
  ! local variables
  INTEGER :: is,lmo,io,jo,i
  INTEGER :: l1,l2,l3,m1,m2,m3
  INTEGER :: lm1,lm2,lm3
  REAL(8) :: t0
  COMPLEX(8) :: z1
  ! allocatable arrays
  COMPLEX(8), ALLOCATABLE :: a(:,:),b(:,:)
=#

    isp = atm2species[ias]
    lmo = lmoapw[isp]
    a = zeros(ComplexF64, lmo, Ngwk)
    b = zeros(ComplexF64, lmo, Ngwk)
    t0 = 0.5d0*rmt[isp]^2
    i = 0
    lm1 = 0
    for l1 in 0:lmaxapw
        for m1 in -l1:l1
            lm1 = lm1 + 1
            for io in 1:apword[isp][l1]
                i = i + 1
                b[i,:] = 0.d0
                lm3 = 0
                for l3 in 0:lmaxapw
                    for m3 in -l3:l3
                        lm3 = lm3 + 1
                        for jo in 1:apword[isp][l3]
                            z1 = 0.0
                            for l2 in 0:lmaxo
                                if mod(l1+l2+l3, 2) == 0 
                                    for m2 in -l2:l2
                                        lm2 = idxlm(l2,m2)
                                        z1 += gntyry(lm2,lm3,lm1)*haa(lm2,jo,l3,io,l1,ias)
                                    end # m2 
                                end 
                            end 
                            if abs(real(z1)) + abs(imag(z1)) > 1.e-14
                                #call zaxpy(ngp, z1, apwalm(:,jo,lm3),1, b(i,1),lmo )
                                b[:,1] .= z1*apwalm[:,jo,lm3] .+ b[:,1]
                            end 
                        end 
                    end 
                end 
                # kinetic surface contribution
                for jo in 1:apword[isp][l1]
                    z1 = t0*apwfr(nrmt(is),1,io,l1,ias)*apwdfr(jo,l1,ias)
                    #CALL zaxpy(ngp,z1,apwalm(:,jo,lm1),1,b(i,1),lmo)
                end 
                a[i,1:Ngwk] = apwalm[1:Ngwk,io,lm1]
            end 
        end 
    end 
    #CALL zmctmu(thr,lmo,ngp,a,b,ld,h)
    return
end
