require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a triangular matrix
NB.
NB. Syntax:
NB.   T=. (uplo ; diag) xtrtri AA
NB. where
NB.   uplo - string, case-insensitive, in which the head
NB.           specifies which triangular part of AA is to be
NB.           referenced:
NB.            'L' - lower
NB.            'U' - upper
NB.   diag - string, case-insensitive, in which the head
NB.           specifies the form of A:
NB.            'N' - A is non-unit triangular: L or U
NB.            'U' - A is unit triangular: L1 or U1
NB.   AA   - n×n-matrix, contains either non-zero or both
NB.          part(s) of A
NB.   T    - n×n-matrix, contains inversed triangular
NB.          matrix in changed triangle and unchanged
NB.          elements of AA in [strict] opposite triangle
NB.   A    - n×n-matrix to be inversed, triangular
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
