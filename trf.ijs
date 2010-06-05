NB. trf.ijs
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
NB. LU factorization of a general M-by-N matrix using partial
NB. pivoting with row interchanges:
NB.   P * L * U = A
NB.
NB. Syntax:
NB.   'p LU'=. getrf A
NB. where
NB.   A  - M-by-N matrix
NB.   LU - M-by-N matrix, the factors L and U from
NB.        factorization
NB.   p  - N-vector, standart permutation vector to pivot
NB.        rows
NB.   L -  M-by-min(M,N) matrix, lower triangular, the unit
NB.        diagonal elements are not stored
NB.   U -  min(M,N)-by-N matrix, upper triangular
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


NB. =========================================================
Note 'trf testing and timing'
   A=: (6 8 ?.@$ 20x) * 0=6 8 ?.@$ 3
   'p LU'=. getrf_pjlap_ A
   L=. ((<./ $ A) ({. " 1) sltri_jlapack_ LU) + idmat_jlapack_ (# A)
   U=. (<./ $ A) {. utri_jlapack_ LU
   P=. p2P_pjlap_ p
   Pinv=. p2P_pjlap_ /: p
   Pinv -: |: P
1
   Pinv -: %. P
1
   A -: ((|: P) mp L mp U)
1
   A -: ((/: p) C. L mp U)
1
   (P mp A) -: (L mp U)
1
   (p C. A) -: (L mp U)
1

   ts=: 6!:2, 7!:2@]
   a500f=. 0.1 * ? 500 500 $ 10
   5 ts 'getrf a500f'
7.13079 9.46784e6
)
