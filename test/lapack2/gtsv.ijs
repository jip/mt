require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a square
NB.   tridiagonal matrix
NB.
NB. Syntax:
NB.   X=. xgtsv dlA ; dA ; duA ; B
NB. where
NB.   dlA  - (n-1)-vector, subdiagonal of A
NB.   dA   - n-vector, diagonal of A
NB.   duA  - (n-1)-vector, superdiagonal of A
NB.   B    - n×nrhs-matrix, RHS
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, will be factored as:
NB.            P * L1 * U = A
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgtsv=: 3 : 0
  'dl d du B'=. y
  'n nrhs'=. $ B
  assert. (isvector_jlapack2_ , isreal_jlapack2_ , (<: n) = #) dl
  assert. (isvector_jlapack2_ , isreal_jlapack2_ ,     n  = #) d
  assert. (isvector_jlapack2_ , isreal_jlapack2_ , (<: n) = #) du
  assert. (ismatrix_jlapack2_ , isreal_jlapack2_             ) B
  cdrc=. dgtsv_jlapack2_ (, n) ; (, nrhs) ; dl ; d ; du ; (|: B) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 6 {:: cdrc
)

zgtsv=: 3 : 0
  'dl d du B'=. y
  'n nrhs'=. $ B
  assert. (isvector_jlapack2_ , (<: n) = #) dl
  assert. (isvector_jlapack2_ ,     n  = #) d
  assert. (isvector_jlapack2_ , (<: n) = #) du
  assert.  ismatrix_jlapack2_               B
  cdrc=. zgtsv_jlapack2_ (, n) ; (, nrhs) ; dl ; d ; du ; (|: B) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 6 {:: cdrc
)
