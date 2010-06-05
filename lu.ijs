NB. lu.ijs
NB. LU factorization
NB.
NB. getrf  LU factorization of a general matrix

coclass 'mt'

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

  mn=. 'm n'=. $ y
  p=. i. m

  for_j. i. m <. n do.
    c2=. j ([ }. ({"1)) y                 NB. col on and under diag
    jp=. j + jp0=. {. \: | c2                   NB. IO max of c2 abs values
NB.smoutput 'getf2fix' ; 'j,jp' ; (j , jp) ; 'p' ; ($p) ; (,.p) ; 'c2' ; ($c2) ; (,.c2) ; 'y' ; ($y) ; y
    if. jp0 do.
      p=. (< jp , j) C. p                NB. accumulate permutations
      y=. (< jp , j) C. y                NB. swap j-th and jp-th rows of y
      c2=. (< 0 , jp0) C. c2             NB. swap j-th and jp-th items of c
NB.smoutput 'jp0≠0' ; 'p' ; ($p) ; p ; 'new c2' ; ($c2) ; c2 ; 'new y' ; ($y) ; y
    end.
    iosr=. m <@(si2i >:) j
    iosc=. n <@(si2i >:) j
    c2=. (}. % {.) c2
NB.c2o=. }. c2
NB.c2=. c2o % {. c2
NB.smoutput 'p' ; ($p) ; (,.p) ; 'old c2' ; ($c2o) ; (,.c2o) ; 'c2' ; ($c2) ; (,.c2) ; 'iosr' ; ($iosr) ; iosr ; 'iosc' ; iosc ; 'y' ; ($y) ; y
    y=. c2 (< iosr , < j) } y
NB.smoutput '((>: j) - mn)' ; (((>: j) - mn)) ; 'y' ; ($y) ; y
    y=. ((((>: j) - mn) {. y) - c2 */ (j (>:@[ }. {) y)) (< iosr , iosc) } y
NB.smoutput 'y' ; ($y) ; y
  end.

  p ; y
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

NB. unblocked version
NB. partial (rows) pivoting
NB. Crout version of the algorithm

NB. NETLIB slus.f

sluc=: (3 : 0) " 2
  si2i=. ] + (i. @ -)                    NB. y (y+1) ... (x-1)

  mn=. <./ 'm n'=. $ y
  p=. i. m

  for_j. i. mn do.
    c=. j ({"1) y                        NB. take j-th col
    c2=. (j }. c) - (((j - m) , j) {. y) mp (j {. c)        NB. compute updated c's tail
    jp=. j + jp0=. {. \: | c2            NB. find pivot
NB.smoutput 'SLUC' ; 'c' ; ($c) ; c ; 'c2' ; ($c2) ; c2 ; 'p' ; ($p) ; p ; 'jp' ; jp
    if. jp0 do.
      p=. (< jp , j) C. p                NB. accumulate permutations
      y=. (< jp , j) C. y                NB. swap j-th and jp-th rows of y
      c2=. (< 0 , jp0) C. c2             NB. swap j-th and jp-th items of c
NB.smoutput 'jp0≠0' ; 'p' ; ($p) ; p ; 'new c2' ; ($c2) ; c2 ; 'new y' ; ($y) ; y
    end.
    c2=. ({. 0 } (% {.)) c2              NB. scale c2's tail by its head
    y=. c2 (< (m si2i j) ; j) } y        NB. write back c2
NB.smoutput 'new c2' ; ($c2) ; c2 ; 'new y' ; ($y) ; y
    r=. j { y
    y=. (((>: j) }. r) - (j {. r) mp (j , (>: j - n)) {. y) (< j ; (n si2i >: j)) } y  NB. write back r2
NB.smoutput 'r' ; ($r) ; r ; 'new y' ; ($y) ; y
  end.

  p ; y
)

slubc=: (3 : 0) " 2
  si2i=. ] + (i. @ -)                    NB. y (y+1) ... (x-1)

  nb=. 32                                NB. choosed experimentally
  mn=. <./ 'm n'=. $ y
  p=. i. m

  for_j. range 0 , (mn - 1) , nb do.
    jb=. nb <. mn - j
    iosc=. (_1 , j) ,: ((m - j) , jb)
    iosa=. (j - m) , j
    iosb=. (0 , j) ,: (j , jb)
    'p2 c'=. sluc ((iosc & (] ;. 0)) - (iosa & {.) mp (iosb & (] ;. 0))) y  NB. update (c←c-a*b), then LU-factorize diagonal and subdiagonal blocks c
    iosr=. m si2i j                      NB. where p2 will modify p
    iosc=. (j + jb) si2i j
    dp=. (i. j) , (j + p2)               NB. current iteration's rows permutation
    y=. dp C. y                          NB. apply rows permutation starting from j-th row
    p=. dp C. p                          NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
    y=. c (< iosr ; iosc) } y            NB. write back diagonal and subdiagonal blocks and hence undo p2 permutation in it
    iosc=. (j , _1) ,: (jb , (n - (j + jb)))
    iosa=. (j , 0) ,: (jb , j)
    iosb=. j , ((j + jb) - n)
    c=. ((iosc & (] ;. 0)) - (iosa & (] ;. 0)) mp (iosb & {.)) y  NB. update (c←c-a*b)
    c=. c getrslx1 ((,.~ j , jb) (] ;. 0) y)                      NB. solve (L*x=c)
    iosr=. (j + jb) si2i j
    iosc=. n si2i (j + jb)
    y=. c (< iosr ; iosc) } y            NB. write back x
  end.

  p ; y
)

NB. Kernel of LU Crout algo
kluc=: 1 : 0
:
  '`f g h i j'=. u
  x ([ (i - (h & (x g y)) @ j) f) y
)

NB. Implement kernel for both directions
kluc_col=: ({"1)`(((- (0 & (1 }) @ $))~ {.)~ {. ])`(mp~)`(( }.   ~ {.)~)`(( {.   ~ {.)~) kluc
kluc_row=:  {   `(((- (0 & (0 }) @ $))~ {.)~ {. ])` mp  `(((}."1)~ {.)~)`((({."1)~ {.)~) kluc

NB. LU Crout algo
NB. 'p LU'=. nb luc A
NB. where nb - block size, nb>0

luc=: (4 : 0) " 0 2
  mkcp=. < @ (, ` (, @ ]) @. =)          NB. make cycle permutation from x and y
  si2i=. ] + (i. @ -)                    NB. y (y+1) ... (x-1)

  mn=. <./ 'm n'=. $ y
  p=. i. m

  if. n > 1 do.
    for_j. range 0 , (mn - 1) , x do.
      jb=. x <. mn - j
      iosc=. (j + jb) si2i j
      'p2 rc2'=. 1 luc iosc kluc_col y   NB. update (c←c-a*b), then factorize diagonal and subdiagonal blocks
      dp=. (i. j) , (j + p2)             NB. current iteration's rows permutation
      p=. dp C. p                        NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
      y=. dp C. y                        NB. apply rows permutation starting from j-th row
      y=. rc2 (< (m si2i j) ; iosc) } y  NB. write back diag. and subdiag. blocks and undo p2 permutation in it
      rc2=. iosc kluc_row y              NB. update (r←r-a*b) diag. and behind diag. blocks
      rc2=. jb ((}. " 1) getrslx1 ({. " 1)) rc2   NB. solve (L*x=c)
      y=. rc2 (< iosc ; (n si2i (j + jb))) } y       NB. write back x
    end.
  else.
    jp=. {. \: | y                       NB. find pivot
    dp=. 0 mkcp jp
    p=. dp C. p                          NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
    y=. ({. 0 } (%"1 {.)) dp C. y        NB. apply rows permutation starting from j-th row, then scale by head
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
  if. (m > 2) *. (n > 2) do. NB. (CPU_CACHE < 7!:5 < 'y') *. (n > 1) do.
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
    sluc y
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
  'p LU'=. sluc y
  smoutput 'sluc (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. 1 luc y
  smoutput '1 & luc (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. slubc y
  smoutput 'slubc (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. 32 luc y
  smoutput '32 & luc (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'L U p'=. getrf_jlapack_ y
  smoutput 'LAPACK getrf (match error): ' , ": y (-: , error) p invperm_jlapack_~ L mp U
NB.  'p L U'=. jLU y
NB.  smoutput 'jwiki LU (match error): ' , ": y (-: , error) p {"1 L mp U
)

NB. measure time and space
NB.   mgetrf mat

mgetrf=: 3 : 0
  ts=. 6!:2, 7!:2@]
  ts & > 'pLU=. getf2 y';'pLU=. getf2fix y';'pLU=. rgetrf y';'pLU=. sluc y';'pLU=. 1 luc y';'pLU=. slubc y';'pLU=. 32 luc y';'LUp=. getrf_jlapack_ y' NB.;'pLU=. jLU y'
)

NB. test and measure interface verbs
NB.   testgetrf m,n
NB.   testgetrf n

testgetrf=: 3 : 0
  y=. 2 {. y , {: y    NB. n->(n,n) , (m,n)->(m,n)
  (tgetrf ] mgetrf) A=: 0.1 * ? y $ 100
)

NB.    testgetrf_mt_ 300 200
NB. getrf (match error): 0 66.3405
NB. getf2r (match error): 0 1.91554e_5
NB. rgetrf (match error): 0 1.91554e_5
NB. LAPACK getrf (match error): 0 1.35102e_5
NB.  1.26323 2.38477e6
NB. 0.622143 1.06835e6
NB. 0.617834 1.06963e6
NB. 0.076975   2.624e6
NB.    testgetrf_mt_ 300 500
NB. |NaN error: getf2r
NB. |   y=.(((i((>:@[)}.({"1))y)-l)    %((n1*i)({,)y))(<(m si2i>:i);i)}y
NB.    testgetrf_mt_ 500 300
NB. |NaN error: getf2r
NB. |   y=.(((i((>:@[)}.({"1))y)-l)    %((n1*i)({,)y))(<(m si2i>:i);i)}y

NB.    testgetrf_mt_ 1000
NB. getf2 (match error): 0 9.69122e_5
NB. getf2fix (match error): 0 9.69122e_5
NB. rgetrf (match error): 0 1459.65
NB. sluc (match error): 0 2.39063e_5
NB. sluc2 (match error): 0 2.39063e_5
NB. slubc (match error): 0 5.42776e_5
NB. LAPACK getrf (match error): 0 9.69768e_5
NB. 58.0344 3.78039e7
NB. 61.1677 3.78041e7
NB. 9.50659 4.19584e7
NB. 41.0749 1.89153e7
NB.  40.391  1.8905e7
NB. 3.19797 1.92858e7
NB. 5.48789 4.08991e7
