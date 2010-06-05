NB. sv.ijs
NB. Solve linear system
NB.
NB. getrsux  Solve system U*X=B, where U is an upper
NB.          triangular matrix
NB. getrslx  Solve system L*X=B, where L is a lower
NB.          triangular matrix
NB. getrsxl  Solve system X*L=B, where L is a lower
NB.          triangular matrix
NB. getrsxu  Solve system X*U=B, where U is an upper
NB.          triangular matrix
NB. getrsux1 Solve system U*X=B, where U is an unit upper
NB.          triangular matrix
NB. getrslx1 Solve system L*X=B, where L is an unit lower
NB.          triangular matrix
NB. getrsxl1 Solve system X*L=B, where L is an unit lower
NB.          triangular matrix
NB. getrsxu1 Solve system X*U=B, where U is an unit upper
NB.          triangular matrix
NB. gesv     Solve system A*X=B via LU factorization with
NB.          partial pivoting
NB. disv     Solve system A*x=b, where A is a diagonalizable
NB.          matrix
NB. hesv     Solve system A*x=b, where A is a Hermitian
NB.          matrix

coclass 'mt'

NB. =========================================================
NB. Local definitions

getrs=: 2 : 0
:
  n1=. >: n=. # y
  for_i. u i. n do.
    x=. (((i { x) - ((i { y) (mp & (i & v)) x)) % ((n1 * i) ({ ,) y)) i } x
  end.
  x
)

getrs1=: 2 : 0
:
  n1=. >: n=. # y
  for_i. u i. n do.
    x=. ((i { x) - ((i { y) (mp & (i & v)) x)) i } x
  end.
  x
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. getrsux                                               2 2
NB. getrslx                                               2 2
NB. getrsxl                                               2 2
NB. getrsxu                                               2 2
NB. Solve systems: U*X=B, L*X=B, X*L=B, X*U=B; where U is an
NB. upper triangular matrix, L is an lower triangular matrix
NB.
NB. Syntax:
NB.   X=. B getrsux LU
NB.   X=. B getrslx LU
NB.   X=. B getrsxl LU
NB.   X=. B getrsxu LU
NB. where
NB.   LU - N×N-matrix, the factors L and U from
NB.        factorization
NB.   B  - N-vector or N×NRHS-matrix, RHS
NB.   X  - same shape as B, solution
NB.   N    >= 0
NB.   NRHS >= 0
NB.
NB. If:<<<<<<<<<<<<<<<<<<
NB.   X1=. B getrsu LU
NB.   U=. utri_jlapack_ LU
NB.   X2=. B getrsu U

NB.   X1=. B getrsl LU
NB.   L=. ltri_jlapack_ LU
NB.   X2=. B getrsl1 L
NB. then
NB.   B -: L mp X1
NB.   X1 -: X2
NB. then
NB.   B -: U mp X1
NB.   X1 -: X2
NB.
NB. Notes:
NB. - strict lower triangle is not referenced
NB. - reuse x to save result
NB. - result is identical to LAPACK's dgetrs/zgetrs
NB. - mp involves excessive pair A[i][i]*z[i] (z[i]==0) to avoid increment
NB.
NB. TODO:
NB. - transform to adverb: ({ getrsu) to process (U*X=B) and ({"1 getrsu) to process (X*L=B)

getrsux=: (|. getrs ((}.~ >:)~)) " 2 2      NB. equivalent to (x mp 128!:1 y)
getrslx=: (] getrs {.) " 2 2
getrsxl=: getrsux &. |:
getrsxu=: getrslx &. |:

NB. ---------------------------------------------------------
NB. getrsux1                                              2 2
NB. getrslx1                                              2 2
NB. getrsxl1                                              2 2
NB. getrsxu1                                              2 2
NB. Solve systems: U*X=B, L*X=B, X*L=B, X*U=B; where U is an
NB. upper triangular matrix with units on diagonal, L is an
NB. lower triangular matrix with units on diagonal
NB.
NB. Syntax:
NB.   X=. B getrsux1 LU
NB.   X=. B getrslx1 LU
NB.   X=. B getrsxl1 LU
NB.   X=. B getrsxu1 LU
NB. where
NB.   LU - N×N-matrix, the factors L and U from
NB.        factorization
NB.   B  - N-vector or N×NRHS-matrix, RHS
NB.   X  - same shape as B, solution
NB.   N    >= 0
NB.   NRHS >= 0
NB.
NB. If:<<<<<<<<<<<<<<<<<<
NB.   X1=. B getrsu LU
NB.   U=. utri_jlapack_ LU
NB.   X2=. B getrsu U
NB. then
NB.   B -: U mp X1
NB.   X1 -: X2

NB.   X1=. B getrsl LU
NB.   L=. ltri_jlapack_ LU
NB.   X2=. B getrsl1 L
NB. then
NB.   B -: L mp X1
NB.   X1 -: X2

NB.   X1=. B getrsl1 LU
NB.   L=. (sltri_jlapack_ LU) + (idmat_jlapack_ # LU)
NB.   X2=. B getrsl1 L
NB. then
NB.   B -: L mp X1
NB.   X1 -: X2

NB.
NB. Notes:
NB. - strict lower triangle is not referenced
NB. - reuse x to save result
NB. - result is identical to LAPACK's dgetrs/zgetrs
NB. - mp involves excessive pair A[i][i]*z[i] (z[i]==0) to avoid increment
NB.
NB. TODO:
NB. - transform to adverb: ({ getrsu) to process (U*X=B) and ({"1 getrsu) to process (X*L=B)

getrsux1=: (|. getrs1 ((}.~ >:)~)) " 2 2
getrslx1=: (] getrs1 {.) " 2 2
getrsxl1=: getrsux1 &. |:
getrsxu1=: getrslx1 &. |:

NB. ---------------------------------------------------------
NB. gesv                                                  1 2
NB. Solve system A*X=B, where A is a general matrix
NB.
NB. Syntax:
NB.   X=. B gesv A
NB. where
NB.   A - N×N-matrix
NB.   B  - N-vector or N×NRHS-matrix, RHS
NB.   X  - same shape as B, solution
NB.   N    >= 0
NB.   NRHS >= 0
NB.
NB. If:
NB.   X=. B gesv A
NB. then
NB.   B -: A mp X
NB.
NB. Notes:
NB. - result is identical to LAPACK's dgesv/zgesv

gesv=: (((0 {:: ]) C. [) (getrslx1 getrsux ]) (1 {:: ])) getrf

NB. ---------------------------------------------------------
NB. disv                                                  1 1
NB. Solve system A*x=b, where A is a diagonalizable matrix
NB.
NB. Syntax:
NB.   x=. b disv (rv ; ev ; rvi)
NB. where
NB.   rv  - N×N-matrix, right eigenvectors of A
NB.   ev  - N-vector, eigenvalues of A
NB.   rvi - N×N-matrix, inversion of rv
NB.   x   - N-vector, solution
NB.   N  >= 0

disv=: [: : ((0 {:: ]) % (1 {:: ])) mp (mp~ (2 & {::)) " 1 1

NB. ---------------------------------------------------------
NB. hesv                                                  1 1
NB. Solve system A*x=b, where A is a Hermitian matrix
NB.
NB. Syntax:
NB.   x=. b hesv (rv ; ev)
NB. where
NB.   rv - N×N-matrix, right eigenvectors of A
NB.   ev - N-vector, eigenvalues of A
NB.   x  - N-vector, solution
NB.   N >= 0

hesv=: [: : ((0 {:: ]) % (1 {:: ])) mp (mp~ (+ @ |: @ (0 & {::))) " 1 1

NB. =========================================================
Note 'trs testing and timing'
   n=. 10
   A=. ? (2 $ n) $ 10
   LU=. 1 {:: getrf_mt_ A
   x=. ? n $ 10
   U=. utri_jlapack_ LU
   bU=. U mp x
   L=. (sltri_jlapack_ LU) + (idmat_jlapack_ n)
   bL=. L mp x
   b=. A mp x
   +/ | x - bU getrsu_mt_ U
1.98175e_14
   +/ | x - bL getrsl1_mt_ L
1.77636e_15
   +/ | x - b gesv_mt_ A
5.86198e_14
   +/ | x - bU trtrsu_mt_ trU
1.98175e_14

   ts=: 6!:2, 7!:2@]
   L1000=. (sltri_jlapack_ 1 {:: getrf_mt_ 0.1 * ? 1000 1000 $ 100) + (idmat_jlapack_ 1000)
   U1000=. utri_jlapack_ 1 {:: getrf_mt_ 0.1 * ? 1000 1000 $ 100
   U1000tr=. ut2tr_mt_ U1000
   x1000=. 0.1 * ? 1000 $ 100
   bU1000=. U1000 mp x1000
   bL1000=. L1000 mp x1000
   10 (ts & >) 'bU1000 getrsu_mt_ U1000';'bL1000 getrsl1_mt_ L1000';'bU1000 trtrsu_mt_ U1000tr';'(128!:1 U1000) mp bU1000'
0.0423727     41984
0.0347479     41600
0.0489226     34176
  1.07046 3.35562e7
   ts 'xx1000=. B1000 gesv_mt_ A1000'
89.0654 3.78045e7
)

NB. =========================================================
NB. Test suite

NB. test for errors
NB.   tgetrs mat

tgetrs=: 3 : 0
  mn=. <./ 'm n'=. $ y  NB. <<<<<<<<<<<<<<<
  im=. i. m
  in=. i. n
  U=. im <:/ in                          NB. upper triangle rectangular matrix
  L=. im >/ in                           NB. strict lower triangle rectangular matrix
  I=. im =/ in                           NB. identity rectangular matrix

  LU2L=. (mn & ({. " 1)) @: (I + L & *)  NB. L is m×min(m,n) matrix
  LU2U=. (mn & {.) @: (U & *)            NB. U is min(m,n)×n matrix
  error=. %: @ (+/) @: , @: | @: -

  'p LU'=. getrf y
  smoutput 'getrf (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. getf2r y
  smoutput 'getf2r (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. rgetrf y
  smoutput 'rgetrf (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
)

NB. measure time and space
NB.   mgetrf mat

mgetrs=: 3 : 0
  ts=. 6!:2, 7!:2@]
  ts & > 'pLU=. getrf y';'pLU=. getf2r y';'pLU=. rgetrf y'
)

NB. test and measure interface verbs
NB.   testgetrf m,n
NB.   testgetrf n

testgetrs=: 3 : 0
  y=. 2 {. y , {: y    NB. n->(n,n) , (m,n)->(m,n)
  (tgetrs ; mgetrs) 0.1 * ? y $ 10
)
