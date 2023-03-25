require 'math/lapack2'

NB. Description:
NB.   Compute the Cholesky factorization of a Hermitian
NB.   (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   T=. uplo xpotrf A
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.            'L' - lower, the form is:
NB.                    L * L^H = A
NB.            'U' - upper, the form is:
NB.                    U^H * U = A
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite or upper or lower triangular
NB.   T    - n×n-matrix, L or U factor combined with
NB.          untouched elements from A in opposite triangle
NB.   L    - n×n-matrix, the lower triangular
NB.   U    - n×n-matrix, the upper triangular
NB.   n    ≥ 0, the size of A, L and U
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - no check for positive definiteness

dpotrf=: 4 : 0
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , isreal_jlapack2_) y
  if. JFL ~: 3!:0 y do. y=. 9 o. y end.
  n=. # y
  cdrc=. dpotrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)

zpotrf=: 4 : 0
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  cdrc=. zpotrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)
