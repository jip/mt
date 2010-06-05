NB. lu.ijs
NB. LU factorization
NB.
NB. getrf  LU factorization of a general matrix

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. getrf                                                   2
NB. LU factorization of a general M×N matrix using partial
NB. pivoting with row interchanges:
NB.   P * L * U = A
NB.
NB. Syntax:
NB.   'p LU'=. getrf A
NB. where
NB.   A  - M×N matrix
NB.   LU - M×N matrix, the factors L and U from
NB.        factorization
NB.   p  - N-vector, standart permutation vector to pivot
NB.        rows
NB.   L -  M×min(M,N) matrix, lower triangular, the unit
NB.        diagonal elements are not stored
NB.   U -  min(M,N)×N matrix, upper triangular
NB.   M >= 0
NB.   N >= 0
NB.
NB. If:
NB.   'p LU'=. getrf A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB.   L=. ((<./ $ A) ({. " 1) sltri_jlapack_ LU) + idmat_jlapack_ (# A)
NB.   U=. (<./ $ A) {. utri_jlapack_ LU
NB. then
NB.   Pinv -: |: P
NB.   Pinv -: p2P /: p
NB.   A -: Pinv mp L mp U
NB.   A -: (/: p) C. L mp U
NB.   (p C. A) -: (L mp U)
NB.
NB. Notes:
NB. - result is identical to LAPACK's dgetrf/zgetrf

NB. LAPACK's getf2

getf2=: (3 : 0) " 2
  mkcp=. < @ (, ` (, @ ]) @. =)          NB. make cycle permutation from x and y
  si2i=. ] + (i. @ -)                    NB. y (y+1) ... (x-1)

  'm n'=. $ y
  p=. iosr=. i. m                        NB. permutations accumulator and row indices
  for_j. i. 0 >. <: m <. n do.           NB. traverse curtailed diagonal, protect against zero-sized y
    jp=. j + {. \: | (< iosr ; j) { y    NB. io max of abs values of diag and under diag element
    cp=. j mkcp jp                       NB. new cycle permutation
    if. ((< jp ; j) { y) ~: 0 do.        NB. if not all elements are zero
      if. j ~: jp do.                    NB. if there is a permutation
        p=. cp C. p                      NB. accumulate permutations
        y=. cp C. y                      NB. pivot row
      end.
      iosr=. }. iosr                     NB. 'rows under diagonal' indices
      ios1=. < iosr ; j                  NB. 'elements in column under diagonal' indices
      sc=. (ios1 { y) % ((< 2 $ j) { y)  NB. scale elements in column under diagonal
      y=. sc ios1 } y                    NB. replace part of column by sc
      iosc=. n si2i >: j                 NB. 'columns behind diagonal' indices
      ios2=. < iosr ; iosc               NB. 'submatrix to adjust' element indices
      y=. ((ios2 { y) - (sc */ ((< j ; iosc) { y))) ios2 } y
    end.
  end.
  p ; y
)

NB. unblocked version
NB. partial (rows) pivoting
NB. right-looking version of the algorithm

NB. LAPACK's getf2 fixed version

getf2fix=: (3 : 0) " 2
  mkcp=. < @ (, ` (, @ ]) @. =)          NB. make cycle permutation from x and y
  si2i=. ] + (i. @ -)                    NB. y (y+1) ... (x-1)

  mn=. <./ 'm n'=. $ y
  n1=. >: n
  p=. iosr=. i. m

  for_i. i. mn do.
    c=. i ([ }. ({"1)) y                 NB. col on and under diag
    ip=. i + {. \: | c                   NB. IO max of c abs values
    cp=. i mkcp ip                       NB. new cycle permutation (no protection against 'all zeros' case)
NB. smoutput 'swap rows: ' ; (i , ip) ; '$p' ; ($p) ; 'c' ; ($c) ; c
    p=. cp C. p                          NB. accumulate permutations
    y=. cp C. y                          NB. swap rows
    iosr=. }. iosr
    ios2=. < iosr ; (n si2i >: i)
    c=. (}. % {.) i ([ }. ({"1)) y
NB. c=. (}. % {.) co=. i ([ }. ({"1)) y
NB.co=. }. co
    y=. c (< iosr ; i) } y
NB.zz=. (ios2 { y)
    y=. ((ios2 { y) - c */ (i (>:@[ }. {) y)) ios2 } y
NB.smoutput 'old c' ; co ; 'c' ; c ; 'A22 before upd' ; ($ zz) ; zz ; 'after' ; ($ ios2 { y) ; (ios2 { y)
  end.

  p ; y 
)

NB. recursive blocked

NB. References:
NB. [1] J. J. Dongarra, S. Hammarling, D. W. Walker. (1996)
NB.     Key Concepts For Parallel Out-Of-Core LU Factorization.
NB.     LAPACK Working Note 110, UT-CS-96-324, April 1996.
NB. [2] J. J. Dongarra, E. F. D'Azevedo. (1997)
NB.     The Design and Implementation of the Parallel
NB.     Out-of-core ScaLAPACK LU, QR, and Cholesky
NB.     Factorization Routines.
NB.     LAPACK Working Note 118, UT-CS-97-347, January 1997.

rgetrf=: (3 : 0) " 2
  'm n'=. mn=. $ y
  if. (CPU_CACHE < 7!:5 < 'y') *. (n > 1) do.
    k=. m <. n2=. >. -: n
    'p1 UL1'=. (0 0 ,: m , k) rgetrf ;. 0 y     NB. factorize A's 1st block column recursively
    UL11=. (0 0 ,: 2 $ k) ] ;. 0 UL1            NB. split A's 1st block column on UL11 and L21
    A21=. (_1 0 ,: (m - k) , k) ] ;. 0 y
    L21=. A21 getrsxu UL11                      NB. solve L21*U11=A21 for L21
    y=. (0 _1 ,: m , (n - k)) (p1 & C.) ;. 0 y  NB. apply p1 to A's 2nd block column
    A12=. (0 0 ,: k , (n - k)) ] ;. 0 y
    U12=. A12 getrslx UL11                      NB. solve L11*U12=A12 for U12
    A22=. (_1 0 ,: (m - k) , (n - k)) ] ;. 0 y
    'p2 UL22'=. rgetrf (A22 - L21 mp U12)       NB. factorize updated A22 recursively
    p=. (k + p2) C. p1                          NB. ((i. k) , (k + p2)) C. p1
    p ; (UL11 ,. U12) , ((p2 C. L21) ,. UL22)
  else.
    getf2fix y
  end.
)

NB. J wiki essay "LU decomposition"

jLU=: 3 : 0
 'm n'=. $ A=. y
 if. 1=m do.
  p ; (=1) ; p{"1 A [ p=. C. (n-1);~.0,(0~:,A)i.1
 else.
  m2=. >.m%2
  'p1 L1 U1'=. jLU m2{.A
  D=. (/:p1) {"1 m2}.A
  F=. m2 {."1 D
  E=. m2 {."1 U1
  FE1=. F mp %. E
  G=. m2}."1 D - FE1 mp U1
  'p2 L2 U2'=. jLU G
  p3=. (i.m2),m2+p2
  H=. (/:p3) {"1 U1
  (p1{p3) ; (L1,FE1,.L2) ; H,(-n){."1 U2
 end.
)

NB. =========================================================
NB. Test suite

NB. test for errors
NB.   tgetrf mat

require '~addons/math/lapack/lapack.ijs'
need_jlapack_ 'getrf'

tgetrf=: 3 : 0
  mn=. <./ 'm n'=. $ y
  im=. i. m
  in=. i. n
  I=. im =/ in                           NB. identity rectangular matrix
  U=. im <:/ in                          NB. upper triangle rectangular matrix
  L=. im >/ in                           NB. strict lower triangle rectangular matrix

  LU2L=. (mn & ({. " 1)) @: (I + L & *)  NB. L is m×min(m,n) matrix
  LU2U=. (mn & {.) @: (U & *)            NB. U is min(m,n)×n matrix
  error=. %: @ (+/) @: , @: | @: -

  'p LU'=. getf2 y
  smoutput 'getf2 (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. getf2fix y
  smoutput 'getf2fix (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. rgetrf y
  smoutput 'rgetrf (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'L U p'=. getrf_jlapack_ y
  smoutput 'LAPACK getrf (match error): ' , ": y (-: , error) p invperm_jlapack_~ L mp U
  'p L U'=. jLU y
  smoutput 'jwiki LU (match error): ' , ": y (-: , error) p {"1 L mp U
)

NB. measure time and space
NB.   mgetrf mat

mgetrf=: 3 : 0
  ts=. 6!:2, 7!:2@]
  ts & > 'pLU=. getf2 y';'pLU=. getf2fix y';'pLU=. rgetrf y';'LUp=. getrf_jlapack_ y';'pLU=. jLU y'
)

NB. test and measure interface verbs
NB.   testgetrf m,n
NB.   testgetrf n

testgetrf=: 3 : 0
  y=. 2 {. y , {: y    NB. n->(n,n) , (m,n)->(m,n)
  (tgetrf ] mgetrf) 0.1 * ? y $ 100
)

NB.    testgetrf_pjlap_ 300 200
NB. getrf (match error): 0 66.3405
NB. getf2r (match error): 0 1.91554e_5
NB. rgetrf (match error): 0 1.91554e_5
NB. LAPACK getrf (match error): 0 1.35102e_5
NB.  1.26323 2.38477e6
NB. 0.622143 1.06835e6
NB. 0.617834 1.06963e6
NB. 0.076975   2.624e6
NB.    testgetrf_pjlap_ 300 500
NB. |NaN error: getf2r
NB. |   y=.(((i((>:@[)}.({"1))y)-l)    %((n1*i)({,)y))(<(m si2i>:i);i)}y
NB.    testgetrf_pjlap_ 500 300
NB. |NaN error: getf2r
NB. |   y=.(((i((>:@[)}.({"1))y)-l)    %((n1*i)({,)y))(<(m si2i>:i);i)}y
