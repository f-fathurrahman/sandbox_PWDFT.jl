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
function wavefmt!(lrstp, ia, atoms, mt_vars, apwlo_vars, ngp, apwalm, evecfv, wfmt)

    atm2species = atoms.atm2species

    lmaxo = mt_vars.lmaxo
    lmaxi = mt_vars.lmaxi
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    nrmti = mt_vars.nrmti
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    npcmt = mt_vars.npcmt
    npcmti = mt_vars.npcmti
    lradstp = mt_vars.lradstp
    idxlm = mt_vars.idxlm

    apword = apwlo_vars.apword
    apwfr = apwlo_vars.apwfr
    nlorb = apwlo_vars.nlorb
    lorbl = apwlo_vars.lorbl
    idxlo = apwlo_vars.idxlo
    lofr = apwlo_vars.lofr

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
    @views wfmt[1:npc] .= 0.0
    # in the original Elk code this array have leading dimension of 2 and of type real(8)
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
                    @views wfmt[1:nrci] .+= z1 * apwfr[ia][l][io][1:nrci,1]
                end
                wfmt[iro:(iro+nrco-1)] .+= z1 * apwfr[ia][l][io][iro:(iro+nrco-1),1]
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
            z1 = evecfv[ngp+idxlo[ia][ilo][lm]]
            if l <= lmaxi
                @views wfmt[1:nrci] .+= z1 * lofr[ia][ilo][1:nrci,1]
            end
            wfmt[iro:(iro+nrco-1)] .+= z1 * lofr[ia][ilo][iro:(iro+nrco-1),1]
        end
    end
    return
end
