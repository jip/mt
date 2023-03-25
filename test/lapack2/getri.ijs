require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a factored general square matrix
NB.
NB. Syntax:
NB.   iA=. xgetri L1U ; ipiv
NB. where
NB.   L1U  - n×n-matrix, L1 and U combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   iA   - n×n-matrix, the inversion of A
NB.   L1   - n×n-matrix, the unit lower triangular (unit
NB.          diagonal not stored)
NB.   U    - n×n-matrix, the upper triangular
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   A    - n×n-matrix, represented in the factored form by
NB.          L1U and ipiv
NB.   n    ≥ 0, the size of A and iA
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgetri=: 3 : 0
  'L1U ipiv'=. y
  n=. # ipiv
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) L1U
  assert.  isvector_jlapack2_                               ipiv
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  cdrc=. dgetri_jlapack2_ (, n) ; (|: L1U) ; (, 1 >. n) ; ipiv ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 2 {:: cdrc
)

zgetri=: 3 : 0
  'L1U ipiv'=. y
  n=. # ipiv
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) L1U
  assert.  isvector_jlapack2_                               ipiv
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  cdrc=. zgetri_jlapack2_ (, n) ; (|: L1U) ; (, 1 >. n) ; ipiv ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 2 {:: cdrc
)
