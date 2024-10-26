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
function wavefmt!(lrstp, ia, ngp, apwalm, evecfv, wfmt)

    isp = atm2species[ia]
    ldi = 2*lmmaxi
    ldo = 2*lmmaxo
    iro = nrmti[isp] + lrstp

    if lrstp == 1
        # not using coarse vs fine grid
        # XXX Are there any use cases for this?
        # XXX Probably for testing
        nrc = nrmt[isp]
        nrci = nrmti[isp]
        npc = npmt[isp]
        npci = npmti[isp]
    elseif lrstp == lradstp
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        npc = npcmt[isp]
        npci = npcmti[isp]
    else
        error("Invalid lrstp = $(lrstp)")
    end
    nrco = nrc - nrci
    #
    # zero the wavefunction
    @views wfmt[:,1:npc] .= 0.0
    #
    # APW functions     
    #
    lm = 0
    for l in 0:lmaxo
        for m in -l:l
            lm = lm + 1
            i = npci + lm
            for io in 1:apword[isp][l]
                z1 = BLAS.dotu(ngp, evecfv, 1, apwalm[:,io,lm], 1)
                if l <= lmaxi
                    @views wfmt[1:nrci,lm] .+= z1 * apwfr[ia][l][io][1:nrci,1]
                end
                wfmt[iro:(iro+nrco-1),i] .+= z1 * apwfr[ia][l][io][iro:(iro+nrco-1),1]
            end
        end
    end
    #
    # local-orbital functions     
    #
    for ilo in 1:nlorb[isp]
        l = lorbl[isp][ilo]
        for m in -l:l
            lm = idxlm[l,m]
            i = npci + lm
            z1 = evecfv[ngp+idxlo[ia][lm,ilo]]
            if l <= lmaxi
                @views wfmt[1:nrci,lm] .+= z1 * lofr[ia][ilo][1:nrci,1]
            end
            wfmt[iro:(iro+nrco-1),i] .+= z1 * lofr[ia][ilo][iro:(iro+nrco-1),1]
        end
    end
    return
end
