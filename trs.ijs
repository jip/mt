NB. trs.ijs
NB. Solve triangular system
NB.
NB. getrsu   Solve system U*x=b, where U is an upper triangular matrix
NB. getrsl1  Solve system L*x=b, where L is a unit lower triangular matrix

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. getrsu                                                1 2
NB. Solve system U*x=b, where U is an upper triangular matrix
NB.
NB. Syntax:
NB.   x=. b getrsu LU
NB. where
NB.   LU - N-by-N matrix, the factors L and U from
NB.        factorization
NB.   b  - N-vector, RHS
NB.   x  - N-vector, solution
NB.   N >= 0
NB.
NB. If:
NB.   x1=. b getrsu LU
NB.   U=. utri_jlapack_ LU
NB.   x2=. b getrsu U
NB. then
NB.   b -: U mp x
NB.   x1 -: x2
NB.
NB. Notes:
NB. - result is identical to LAPACK's dgetrs/zgetrs
NB. - mp involves excessive pair A[i][i]*z[i] (z[i]==0) to avoid increment

getrsu=: (4 : 0) " 1 2
  n=. # y
  z=. n $ 0
  for_i. |. i. n do.           NB. traverse all rows of U
    z=. (((i { x) - ((i }. i { y) mp (i }. z))) % ((< 2 $ i) { y)) i } z
  end.
)

NB. ---------------------------------------------------------
NB. getrsl1                                               1 2
NB. Solve system L*x=b, where L is a lower triangular matrix
NB. with units on diagonal
NB.
NB. Syntax:
NB.   x=. b getrl1 LU
NB. where
NB.   LU - N-by-N matrix, the factors L and U from
NB.        factorization
NB.   b  - N-vector, RHS
NB.   x  - N-vector, solution
NB.   N >= 0
NB.
NB. If:
NB.   x1=. b getrsl1 LU
NB.   L=. (sltri_jlapack_ LU) + (idmat_jlapack_ # LU)
NB.   x2=. b getrsl1 L
NB. then
NB.   b -: L mp x
NB.   x1 -: x2
NB.
NB. Notes:
NB. - result is identical to LAPACK's dgetrs/zgetrs

getrsl1=: (4 : 0) " 1 2
  n=. # y
  z=. n $ 0
  for_i. i. n do.              NB. traverse all rows of L
    z=. ((i { x) - ((i {. i { y) mp (i {. z))) i } z
  end.
)

NB. =========================================================
Note 'trs testing and timing'
   n=. 10
   A=. ? (2 $ n) $ 10
   LU=. 1 {:: getrf_pjlap_ A
   x=. ? n $ 10
   U=. utri_jlapack_ LU
   bU=. U mp x
   L=. (sltri_jlapack_ LU) + (idmat_jlapack_ n)
   bL=. L mp x
   x -: bU getrsu_pjlap_ U
1
   x -: bL getrsl1_pjlap_ L
1

   L1000=. (sltri_jlapack_ 1 {:: getrf_pjlap_ 0.1 * ? 1000 1000 $ 100) + (idmat_jlapack_ 1000)
   U1000=. utri_jlapack_ 1 {:: getrf_pjlap_ 0.1 * ? 1000 1000 $ 100
   x1000=. 0.1 * ? 1000 $ 100
   b1000=. U1000 mp x1000
   10 ts 'b1000 getrsu_pjlap_ U1000'
0.0334183 41984
   10 ts 'b1000 getrsl1_pjlap_ L1000'
0.0268176 41344
)
