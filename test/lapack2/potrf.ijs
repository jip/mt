require 'math/lapack2'

NB. Description:
NB.   Compute the Cholesky factorization of a symmetric
NB.   (Hermitian) positive definite matrix
NB.
NB. Syntax:
NB.   T=. uplo xpotrf A
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - to factorize A as: L * L^H = A, and
NB.                  strict lower triangle of A is not
NB.                  referenced
NB.            'U' - to factorize A as: U^H * U = A, and
NB.                  strict upper triangle of A is not
NB.                  referenced
NB.   A    - n×n-matrix, contains triangular matrix to
NB.          factorize
NB.   T    - n×n-matrix, L or U factor
NB.   L    - n×n-matrix, lower triangular
NB.   U    - n×n-matrix, upper triangular
NB.   n    ≥ 0, size of A, L and U
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - no check for positive definiteness

dpotrf=: 4 : 0
  n=. # y
  assert. *./ (1 = # x ) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. issymmetric_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  cdrc=. dpotrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)

zpotrf=: 4 : 0
  n=. # y
  assert. *./ (1 = # x ) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. ishermitian_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  cdrc=. zpotrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)
