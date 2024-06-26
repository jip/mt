require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a square
NB.   matrix using the factorization computed by xGETRF
NB.
NB. Syntax:
NB.   X=. trans xgetrs L1U ; ipiv ; B
NB. where
NB.   trans - string, case-insensitive, in which the head
NB.           specifies the form of the system of equations:
NB.             'N' - A   * X = B  (no transpose)
NB.             'T' - A^T * X = B  (transpose)
NB.             'C' - A^H * X = B  (conjugate transpose)
NB.   L1U   - n×n-matrix, L1 and U combined
NB.   ipiv  - n-vector, integer, pivot indices that define P
NB.   B     - n×nrhs-matrix, RHS
NB.   X     - n×nrhs-matrix, solutions of equation:
NB.             op(A) * X = B
NB.   A     - n×n-matrix, represented in the factored form by
NB.           L1U and ipiv
NB.   L1    - n×n-matrix, unit lower triangular (unit
NB.           diagonal is not stored)
NB.   U     - n×n-matrix, upper triangular
NB.   P     - n×n-matrix, boolean, the permutation matrix
NB.   n     ≥ 0, the order of system
NB.   nrhs  ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgetrs=: 4 : 0
  'L1U ipiv B'=. y
  'n nrhs'=. $ B
  assert. 'nNtTcC' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) L1U
  assert. (isvector_jlapack2_ ,                      n = #) ipiv
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  |: 7 {:: dgetrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: L1U) ; ld ; ipiv ; (|: B) ; ld ; , _1
)

zgetrs=: 4 : 0
  'L1U ipiv B'=. y
  'n nrhs'=. $ B
  assert. 'nNtTcC' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) L1U
  assert. (isvector_jlapack2_ ,                      n = #) ipiv
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  |: 7 {:: zgetrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: L1U) ; ld ; ipiv ; (|: B) ; ld ; , _1
)
