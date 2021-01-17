require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a triangular matrix
NB.
NB. Syntax:
NB.   T=. uplo xtrtri A
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only
NB.            'U' - use upper triangle of A only
NB.   diag - scalar, character, case-insensitive:
NB.            'N' - A is non-unit triangular
NB.            'U' - A is unit triangular, the diagonal
NB.                  elements of A are not
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
  assert. *./ (1 = # uplo) , uplo  e. 'lLuU'
  assert. *./ (1 = # diag) , diag  e. 'nNuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  n=. # y
  cdrc=. dtrtri_jlapack2_ (, uplo) ; (, diag) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 4 {:: cdrc
)

ztrtri=: 4 : 0
  'uplo diag'=. x
  assert. *./ (1 = # uplo) , uplo  e. 'lLuU'
  assert. *./ (1 = # diag) , diag  e. 'nNuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  cdrc=. ztrtri_jlapack2_ (, uplo) ; (, diag) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 4 {:: cdrc
)
