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

getrf=: (3 : 0) " 2
  mkcp=. < @ (, ` (, @ ]) @. =)          NB. make cycle permutation from x and y
  si2i=. ] + (i. @ -)                    NB. y (y+1) ... (x-y-1)

  'm n'=. $ y
  p=. iosr=. i. m                        NB. permutations accumulator and row indices
  for_j. i. 0 >. <: m <. n do.           NB. traverse curtailed diagonal, protect against zero-sized y
    jp=. j + {. \: | (< iosr ; j) { y    NB. io max of abs values of diag and under diag element
    cp=. j mkcp jp                       NB. new cycle permutation
    if. ((< jp ; j) { y) ~: 0 do.        NB. if not all elements are zero
      if. j ~: jp do.                    NB. if there is a permutation
NB.smoutput 'swap rows: ' , ": j , jp
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

NB. unblocked
NB. partial (rows) pivoting
getf2r=: (3 : 0) " 2
  mkcp=. < @ (, ` (, @ ]) @. =)          NB. make cycle permutation from x and y
  si2i=. ] + (i. @ -)                    NB. y (y+1) ... (x-y-1)

  mn=. <./ 'm n'=. $ y
  n1=. >: n
  p=. i. mn

  for_i. i. mn do.
    c=. i ([ }. ({"1)) y                 NB. col on and under diag
    ip=. i + {. \: | c                   NB. IO max of c abs values
    cp=. i mkcp ip                       NB. new cycle permutation (no protection against 'all zeros' case)
    if. i ~: ip do.                      NB. if there is a permutation
smoutput 'swap rows: ' ; (i , ip) ; '$p' ; $p
      p=. cp C. p                        NB. accumulate permutations
      y=. cp C. y                        NB. swap rows (fixme: touch A22 only!)
    end.
    u=. i (([ }. {) - (([ {. {) mp ((0 _1 ,: [ , n - [) (] ;. 0) ]))) y
    y=. u (< i ; n si2i i) } y
NB.smoutput 'i' ; i ; 'u' ; ($u) ; u ; 'c' ; ($c) ; c ; 'cp' ; (<cp) ; 'y' ; ($y) ; y
    l=. ((_1 0 ,: (<: m - i) , i) ] ;. 0 y) mp (i ([ {. ({"1)) y)
    y=. (((i ((>:@[) }. ({"1)) y) - l) % ((n1 * i) ({ ,) y)) (< (m si2i >: i) ; i) } y
NB.smoutput 'l' ; ($l) ; l ; 'y' ; ($y) ; y
  end.

  p ; y 
)

NB. recursive blocked
NB. S. Toledo, Locality of reference in lu decomposition with partial pivoting, Tech. Report RC
NB. 20344(1/19/96), IBM Research Division, T. J. Watson Research Center, Yorktown Heights, New York, 1996.

rgetrf=: (3 : 0) " 2
  if. CPU_CACHE < 7!:5 < 'y' do.
    'm2 n2'=. m2n2=. >. -: 'm n'=. mn=. $ y
    p2q2=. 'p2 q2'=. mn-m2n2
    'p1 UL11'=. (0 0 ,: m2n2) rgetrf ;. 0 y                         NB. factorize A11 recursively
    UL21=. ((_1 0 ,: p2 , n2) ] ;. 0 y) getrsxu UL11                NB. solve L21*U11=A21 for L21
    y=. (0 _1 ,: m , q2) (p1 & C.) ;. 0 y                           NB. apply p1 to A's 2nd block column
    UL12=. ((0 0 ,: m2 , q2) ] ;. 0 y) getrslx UL11                 NB. solve L11*U12=A12 for permuted U12
    'p2 UL22'=. rgetrf (((_1 0 ,: p2q2) ;. 0 y) - UL21 mp UL12)  NB. factorize updated A22 recursively
    (p2 C. p1) ; (UL11 ,. UL12) , (UL21 ,. UL22)
  else.
    getf2r y
  end.
)

NB. =========================================================
NB. testing and timing

testgetrf=: 3 : 0
  mn=. <./ 'm n'=. $ y
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

timegetrf=: 3 : 0
  ts=. 6!:2, 7!:2@]
  ts & > 'pLU=. getrf y';'pLU=. getf2r y';'pLU=. rgetrf y'
)

testtimegetrf=: 3 : 0
  y=. 2 {. y , {: y    NB. n->(n,n) , (m,n)->(m,n)
  a=. 0.1 * ? y $ 10
  testgetrf a
  timegetrf a
)
