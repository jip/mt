require 'math/lapack2'

NB. Description:
NB.   Balance a pair of general square matrices
NB.
NB. Syntax:
NB.   'Abal Bbal ilo ihi lscale rscale'=. job xggbal A ;  B
NB.   'Abal Bbal ilo ihi lscale rscale'=. job xggbal A ,: B
NB. where
NB.   job       - literal, case-insensitive, in which the
NB.               head specifies what to do:
NB.                 'N' - to init ilo, ihi and xscale
NB.                 'P' - to permute only
NB.                 'S' - to scale only
NB.                 'B' - to do both permute and scale
NB.   A,B       - n×n-matrix, a matrix pair
NB.   Abal,Bbal - n×n-matrix, a balanced matrix pair (A,B)
NB.   ilo       ∈ [1,max(1,ihi)], IO starting row and column,
NB.               1-based
NB.   ihi       ∈ [min(ilo,n),n], IO ending row and column,
NB.               1-based
NB.   lscale    - n-vector, float, details of the
NB.               permutations and scaling factors applied
NB.               to the left side of A and B
NB.   rscale    - n-vector, float, details of the
NB.               permutations and scaling factors applied
NB.               to the right side of A and B
NB.   n         ≥ 0, the size of A, B, Abal, Bbal, lscale and
NB.               rscale
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dggbal=: 4 : 0
  'A B'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_        ) A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) B
  ld=. , 1 >. n
  work=. ((6 * n) >.^:('sSbB' e.~ {. x) 1) $ 0.0
  'A B ilo ihi lscale rscale'=. 3 5 7 8 9 10 { dggbal_jlapack2_ (, x) ; (, n) ; (|: A) ; ld ; (|: B) ; ld ; (, 0) ; (, 0) ; (n $ 0.0) ; (n $ 0.0) ; work ; , _1
  (|: A) ; (|: B) ; ({. ilo) ; ({. ihi) ; lscale ; rscale
)

zggbal=: 4 : 0
  'A B'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_        ) A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) B
  ld=. , 1 >. n
  work=. ((6 * n) >.^:('sSbB' e.~ {. x) 1) $ 0.0
  'A B ilo ihi lscale rscale'=. 3 5 7 8 9 10 { zggbal_jlapack2_ (, x) ; (, n) ; (|: A) ; ld ; (|: B) ; ld ; (, 0) ; (, 0) ; (n $ 0.0) ; (n $ 0.0) ; work ; , _1
  (|: A) ; (|: B) ; ({. ilo) ; ({. ihi) ; lscale ; rscale
)
