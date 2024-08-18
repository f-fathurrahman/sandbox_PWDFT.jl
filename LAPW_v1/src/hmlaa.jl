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
function hmlaa!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, haa, H)
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

    atm2species = atoms.atm2species
    
    Ngwk = pw.gvecw.Ngw[ik]
    
    rmt = mt_vars.rmt
    lmaxapw = mt_vars.lmaxapw
    lmaxo = mt_vars.lmaxo
    idxlm = mt_vars.idxlm
    gntyry = mt_vars.gntyry
    nrmt = mt_vars.nrmt

    lmoapw = apwlo_vars.lmoapw
    apword = apwlo_vars.apword
    apwfr = apwlo_vars.apwfr
    apwdfr = apwlo_vars.apwdfr

    isp = atm2species[ia]
    lmo = lmoapw[isp]
    a = zeros(ComplexF64, lmo, Ngwk)
    b = zeros(ComplexF64, lmo, Ngwk)
    t0 = 0.5*rmt[isp]^2
    @info "t0 = $t0"
    i = 0
    lm1 = 0
    for l1 in 0:lmaxapw, m1 in -l1:l1
        lm1 += 1
        for io in 1:apword[isp][l1]
            i += 1
            b[i,:] .= 0.0
            lm3 = 0
            for l3 in 0:lmaxapw, m3 in -l3:l3
                lm3 += 1
                for jo in 1:apword[isp][l3]
                    z1 = 0.0 + im*0.0 # complex
                    for l2 in 0:lmaxo
                        if mod(l1+l2+l3, 2) == 0 
                            for m2 in -l2:l2
                                lm2 = idxlm[l2,m2]
                                z1 += gntyry[lm2,lm3,lm1]*haa[ia][lm2,jo,l3,io,l1]
                            end # m2 
                        end 
                    end 
                    #println("z1 = ", z1)
                    if abs(real(z1)) + abs(imag(z1)) > 1.e-14
                        #call zaxpy(ngp, z1, apwalm(:,jo,lm3),1, b(i,1),lmo )
                        b[i,1:Ngwk] .+= z1*apwalm[ia][1:Ngwk,jo,lm3]
                    end # if
                end # jo
            end # l3, m3
            #
            # kinetic surface contribution
            for jo in 1:apword[isp][l1]
                z1 = t0*apwfr[ia][l1][io][nrmt[isp],1] * apwdfr[ia][l1][jo]
                #println("z1 = ", z1)
                #CALL zaxpy(ngp,z1,apwalm(:,jo,lm1),1,b(i,1),lmo)
                b[i,1:Ngwk] .+= z1*apwalm[ia][1:Ngwk,jo,lm1]
            end 
            a[i,1:Ngwk] = apwalm[ia][1:Ngwk,io,lm1]
        end  # io
    end  # l1, m1
    println("Before zmctmu: sum(a) = ", sum(a))
    println("Before zmctmu: sum(b) = ", sum(b))
    zmctmu!(a, b, H)
    return
end
