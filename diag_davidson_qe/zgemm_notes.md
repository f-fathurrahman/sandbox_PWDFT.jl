`
C := alpha*op( A )*op( B ) + beta*C,
``

(M x N) = (M x K ) (K x N)

[in]    M   

          M is INTEGER
           On entry,  M specifies the number of rows  of the  matrix
           op( A )  and of the  matrix  C.  M must be at least  zero.

[in]    N   

          N is INTEGER
           On entry,  N  specifies the number  of columns of the matrix
           op( B ) and the number of columns of the matrix C. N must be
           at least zero.

[in]    K   
          K is INTEGER
           On entry,  K  specifies  the number of columns of the matrix
           op( A ) and the number of rows of the matrix op( B ). K must
           be at least  zero.
