function zmctmu!(tcr, a, b, c)
#=
  IMPLICIT NONE 
  ! arguments
  logical, intent(in) :: tcr
  INTEGER, intent(in) :: l,n
  COMPLEX(8), intent(in) :: a(l,n),b(l,n)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(out) :: c(*)
  ! local variables
  INTEGER :: l2,i,j,k
  ! external functions
  REAL(8) ddot
  COMPLEX(8) zdotc
  external ddot,zdotc
=#
    l = size(a, 1)
    n = size(a, 2)
    @assert l == size(b, 1)
    @assert n = size(b, 2)
    ld = size(c, 1)

    # XXX need tcr? or use dot for both case

    if tcr 
        # matrix c is real valued
        #l2 = 2*l  # should use 2*l
        for j in 1:n
            k = (j-1)*ld
            for i in 1:j
                k = k + 1
                c[k] += dot( a[:,i], b[:,j] )
            end 
        end 
    else 
        # matrix c is complex valued
        # This is the usual case
        for j in 1:n
            k = (j-1)*ld
            for i in 1:j
                k = k + 1
                c[k] += dot( a[:,i], b[:,j] )
            end
        end
    end
    return
end
