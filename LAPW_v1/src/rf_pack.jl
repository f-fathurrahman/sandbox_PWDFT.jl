function rf_pack(tpack,n,np,ld,rfmt,rfir,v)
#=
  use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: np(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rfmt(ld,natmtot),rfir(ngtot)
real(8), intent(out) :: v(*)
! local variables
integer is,ias
=#

    # pack the function
    for ia in 1:Natoms
        isp = atm2species[ia]
        #dcopy( np(is), rfmt(:,ias),1, v(n+1),1)
        n = n + np[is]
    end
    #call dcopy(ngtot, rfir,1, v(n+1), 1)
    n += ngtot
    return n
return
end


