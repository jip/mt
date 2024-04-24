require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a square
NB.   matrix, and factor the last one
NB.
NB. Syntax:
NB.   'L1U ipiv X'=. xgesv A ; B
NB. where
NB.   A    - n×n-matrix to be factored to L1U and ipiv as:
NB.            P * L1 * U = A
NB.   B    - n×nrhs-matrix, RHS
NB.   L1U  - n×n-matrix, L1 and U combined
NB.   ipiv - n-vector, integer, pivot indices that define P
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   L1   - n×n-matrix, unit lower triangular (unit diagonal
NB.          is not stored)
NB.   U    - n×n-matrix, upper triangular
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgesv=: 3 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) A
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  (|:L:0) 3 5 6 { dgesv_jlapack2_ (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; , _1
    NB. (|:) doesn't affect to ipiv
)

zgesv=: 3 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) A
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  (|:L:0) 3 5 6 { zgesv_jlapack2_ (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; , _1
    NB. (|:) doesn't affect to ipiv
)
