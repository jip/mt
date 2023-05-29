require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a triangular matrix
NB.
NB. Syntax:
NB.   T=. (uplo ; diag) xtrtri A
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.            'L' - lower
NB.            'U' - upper
NB.   diag - literal, case-insensitive, in which the head
NB.           specifies the form of A:
NB.            'N' - A is non-unit triangular
NB.            'U' - A is unit triangular, diagonal
NB.                  elements of A are not referenced
NB.   A    - n×n-matrix to be inversed, [unit]
NB.          {lower,upper}-triangular or general which
NB.          contains the triangle to invert
NB.   T    - n×n-matrix, contains inversed triangular
NB.          matrix in changed triangle and unchanged
NB.          elements of A in [strict] opposite triangle
NB.   n    ≥ 0, the size of A and T
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dtrtri=: 4 : 0
  'uplo diag'=. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  |: 4 {:: dtrtri_jlapack2_ (, uplo) ; (, diag) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
)

ztrtri=: 4 : 0
  'uplo diag'=. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  |: 4 {:: ztrtri_jlapack2_ (, uplo) ; (, diag) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
)
