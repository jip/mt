require 'math/lapack2'

NB. Description:
NB.   Balance a square matrix
NB.
NB. Syntax:
NB.   'Abal ilo ihi scale'=. job xgebal A
NB. where
NB.   job   - string, case-insensitive, in which the head
NB.           specifies what to do:
NB.             'N' - to init ilo, ihi and scale
NB.             'P' - to permute only
NB.             'S' - to scale only
NB.             'B' - to do both permute and scale
NB.   A     - n×n-matrix to balance
NB.   Abal  - n×n-matrix, a balanced matrix A
NB.   ilo   ∈ [1,max(1,ihi)], IO starting row and column,
NB.           1-based
NB.   ihi   ∈ [min(ilo,n),n], IO ending row and column,
NB.           1-based
NB.   scale - n-vector, real, details of the permutations
NB.           and scaling factors applied to A
NB.   n     ≥ 0, the size of A, Abal and scale
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgebal=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  'y ilo ihi scale'=. 3 5 6 7 { dgebal_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (, 0) ; (, 0) ; (n $ 0.0) ; , _1
  (|: y) ; ({. ilo) ; ({. ihi) ; scale
)

zgebal=: 4 : 0
  assert. 'nNpPsSbB' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  'y ilo ihi scale'=. 3 5 6 7 { zgebal_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (, 0) ; (, 0) ; (n $ 0.0) ; , _1
  (|: y) ; ({. ilo) ; ({. ihi) ; scale
)
