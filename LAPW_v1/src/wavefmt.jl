#=
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients (in,complex(ngkmax,apwordmax,lmmaxapw))
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   wfmt   : complex muffin-tin wavefunction passed in as real array
!            (out,real(2,*))

integer, intent(in) :: lrstp,ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax)
real(8), intent(out) :: wfmt(2,*)

! local variables
integer is,ldi,ldo,io,ilo
integer nrc,nrci,nrco,iro
integer l,m,lm,npc,npci,i
complex(8) z1
=#
function wavefmt(lrstp, ia, ngp, apwalm, evecfv, wfmt)

    isp = atm2species[ia]
    ldi = 2*lmmaxi
    ldo = 2*lmmaxo
    iro = nrmti[isp] + lrstp

    if lrstp == 1
        nrc = nrmt(is)
        nrci = nrmti(is)
        npc = npmt(is)
        npci = npmti(is)
    elseif lrstp == lradstp
        nrc = nrcmt(is)
        nrci = nrcmti(is)
  npc=npcmt(is)
  npci=npcmti(is)
else
  write(*,*)
  write(*,'("Error(wavefmt): invalid lrstp : ",I8)') lrstp
  write(*,*)
  stop
end if
nrco = nrc - nrci

! zero the wavefunction
wfmt(:,1:npc)=0.d0

!-----------------------!
!     APW functions     !
!-----------------------!
lm=0
do l=0,lmaxo
  do m=-l,l
    lm=lm+1
    i=npci+lm
    do io=1,apword(l,is)
      z1=zdotu(ngp,evecfv,1,apwalm(:,io,lm),1)
      if (abs(dble(z1)) > 1.d-14) then
        if (l <= lmaxi) then
          call daxpy(nrci,dble(z1),apwfr(1,1,io,l,ias),lrstp,wfmt(1,lm),ldi)
        end if
        call daxpy(nrco,dble(z1),apwfr(iro,1,io,l,ias),lrstp,wfmt(1,i),ldo)
      end if
      if (abs(aimag(z1)) > 1.d-14) then
        if (l <= lmaxi) then
          call daxpy(nrci,aimag(z1),apwfr(1,1,io,l,ias),lrstp,wfmt(2,lm),ldi)
        end if
        call daxpy(nrco,aimag(z1),apwfr(iro,1,io,l,ias),lrstp,wfmt(2,i),ldo)
      end if
    end do
  end do
end do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    i=npci+lm
    z1=evecfv(ngp+idxlo(lm,ilo,ias))
    if (abs(dble(z1)) > 1.d-14) then
      if (l <= lmaxi) then
        call daxpy(nrci,dble(z1),lofr(1,1,ilo,ias),lrstp,wfmt(1,lm),ldi)
      end if
      call daxpy(nrco,dble(z1),lofr(iro,1,ilo,ias),lrstp,wfmt(1,i),ldo)
    end if
    if (abs(aimag(z1)) > 1.d-14) then
      if (l <= lmaxi) then
        call daxpy(nrci,aimag(z1),lofr(1,1,ilo,ias),lrstp,wfmt(2,lm),ldi)
      end if
      call daxpy(nrco,aimag(z1),lofr(iro,1,ilo,ias),lrstp,wfmt(2,i),ldo)
    end if
  end do
end do
return
end subroutine
