function calc_match_coeffs!(ngp, vgpc, gpc, sfacgp, apwalm)
    #=
  ! arguments
  INTEGER, intent(in) :: ngp
  REAL(8), intent(in) :: vgpc(3,ngkmax),gpc(ngkmax)
  COMPLEX(8), intent(in) :: sfacgp(ngkmax,natmtot)
  COMPLEX(8), intent(out) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  ! local variables
  INTEGER :: is,ia,ias,omax
  INTEGER :: l,m,lm,io,jo,i
  INTEGER :: nr,ir,igp,info
  REAL(8) :: t0,t1
  COMPLEX(8) :: z1,z2,z3
  ! automatic arrays
  INTEGER :: ipiv(apwordmax)
  COMPLEX(8) :: a(apwordmax,apwordmax)
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: djl(:,:,:)
  COMPLEX(8), ALLOCATABLE :: ylmgp(:,:),b(:,:)
  ! external functions
  REAL(8) polynm
  external polynm
=#
    djl = OffsetArray(
        zeros(Float64, lmaxapw+1, apwordmax, ngp),
        0:lmaxapw, 1:apwordmax, 1:ngp
    )
    ylmgp = zeros(ComplexF6, lmmaxapw, ngp)

    if apwordmax > 1
        b = zeros(ComplexF64, apwordmax, ngp*(2*lmaxapw+1))
    end
  
    # compute the spherical harmonics of the G+p-vectors
    for igp in 1:ngp
        @views genylmv!( lmaxapw, vgpc[:,igp], ylmgp[:,igp] )
    end

    t0 = 4π/sqrt(CellVolume)
  
    # loop over species
    for isp in 1:Nspecies
        nr = nrmt(is)
        # maximum APW order for this species
        omax = maximsum( specs_info[isp].apword )
        # special case of omax=1
        if omax==1 
            for igp in 1:ngp
                t1 = gpc[igp]*rmt[isp]
                @views sbessel!(lmaxapw, t1, djl[:,1,igp] )
            end
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                for l in 0:lmaxapw
                    z1=(t0/apwfr(nr,1,1,l,ias))*zil(l)
                    for igp in 1:ngp
                        z2 = djl[l,1,igp] * z1 * sfacgp[igp,ia]
                        for m in -l:l
                            lm = idxlm(l,m)
                            apwalm[ia][igp,1,lm] = z2*conjg(ylmgp[lm,igp])
                        end 
                    end 
                end 
            end 
            continue # CYCLE  # next species
        end # if omax == 1 
        # starting point on radial mesh for fitting polynomial of order npapw
        ir = nr-npapw+1
        # evaluate the spherical Bessel function derivatives for all G+p-vectors
        for igp in 1:ngp
            t1 = gpc(igp)*rmt(is)
            for io in 1:omax
                @views sbesseldm!( io-1, lmaxapw, t1, djl[:,io,igp] )
            end 
            t1 = 1.0
            for io in 2:omax
                t1 = t1*gpc[igp]
                djl[:,io,igp] = t1*djl[:,io,igp]
            end 
        end 
        # loop over atoms
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            # begin loop over l
            for l in 0:lmaxapw
                z1 = t0*im^l  #zil[l]
                # set up matrix of derivatives
                for jo in 1:apword[isp][l], io in 1:apword[isp][l]
                    a(io,jo)=polynm(io-1,npapw,rsp(ir,is),apwfr(ir,1,jo,l,ias),rmt(is))
                end 
                # set up target vectors
                i = 0
                for igp in 1:ngp
                    z2=z1*sfacgp(igp,ias)
                    for m in -l:l
                        lm = idxlm(l,m)
                        i = i + 1
                        z3 = z2*conjg(ylmgp[lm,igp])
                        for io in 1:apword[isp][l]
                            b[io,i] = djl[l,io,igp]*z3
                        end 
                    end 
                end 
                # solve the general complex linear systems
                #CALL zgesv(apword(l,is),i,a,apwordmax,ipiv,b,apwordmax,info)
                #
                i = 0
                for igp in 1:ngp
                    for m in -l:l
                        lm = idxlm[l,m]
                        i = i + 1
                        for io in 1:apword[isp][l]
                            apwalm[ia][igp,io,lm] = b[io,i]
                        end
                    end 
                end 
            end  # end loop over l 
        end 
    end # end loops over atoms and species

    return
end

