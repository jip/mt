require 'math/lapack2'

NB. Description:
NB.   Solve a triangular system:
NB.     op(A) * X = B
NB.
NB. Syntax:
NB.   X=. (uplo ; trans ; diag) xtrtrs AA ; B
NB. where
NB.   uplo  - literal, case-insensitive, in which the head
NB.            specifies which triangular part of AA is to be
NB.            referenced:
NB.             'L' - lower
NB.             'U' - upper
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N' - op(A) := A    (no transpose)
NB.             'T' - op(A) := A^T  (transpose)
NB.             'C' - op(A) := A^H  (conjugate transpose)
NB.   diag  - literal, case-insensitive, in which the head
NB.            specifies the form of A:
NB.             'N' - A is non-unit triangular: L or U
NB.             'U' - A is unit triangular: L1 or U1
NB.   AA    - n×n-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   A     - n×n-matrix, triangular
NB.   B     - n×nrhs-matrix, RHS
NB.   X     - n×nrhs-matrix, solutions
NB.   n     ≥ 0, the order of system
NB.   nrhs  ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale
NB. - trans='T' and trans='C' are identic for dtrtrs

dtrtrs=: 4 : 0
  'uplo trans diag'=. x
  'AA B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) AA
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  |: 8 {:: dtrtrs_jlapack2_ (, uplo) ; (, trans) ; (, diag) ; (, n) ; (, nrhs) ; (|: AA) ; ld ; (|: B) ; ld ; , _1
)

ztrtrs=: 4 : 0
  'uplo trans diag'=. x
  'AA B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) AA
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  |: 8 {:: ztrtrs_jlapack2_ (, uplo) ; (, trans) ; (, diag) ; (, n) ; (, nrhs) ; (|: AA) ; ld ; (|: B) ; ld ; , _1
)
