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
NB.            'U' - A is unit triangular, the diagonal
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
NB. - the verbs below are loaded into the current locale

dtrtri=: 4 : 0
  'uplo diag'=. x
  assert. 'lLuU' e.~ {. uplo
  assert. 'nNuU' e.~ {. diag
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , isreal_jlapack2_) y
  if. JFL ~: 3!:0 y do. y=. 9 o. y end.
  n=. # y
  cdrc=. dtrtri_jlapack2_ (, uplo) ; (, diag) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 4 {:: cdrc
)

ztrtri=: 4 : 0
  'uplo diag'=. x
  assert. 'lLuU' e.~ {. uplo
  assert. 'nNuU' e.~ {. diag
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  cdrc=. ztrtri_jlapack2_ (, uplo) ; (, diag) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 4 {:: cdrc
)
