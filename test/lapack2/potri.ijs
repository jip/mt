require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a factored Hermitian (symmetric)
NB.   positive definite matrix
NB.
NB. Syntax:
NB.   iAA=. uplo xpotri T
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    L * L^H = A
NB.            'U' - use upper triangle of A only, form is:
NB.                    U^H * U = A
NB.   T    - n×n-matrix, contains the triangle factor either
NB.          L or U to invert
NB.   iAA  - n×n-matrix, contains the triangular part of iA
NB.          in changed triangle and unchanged elements of T
NB.          in opposite strict triangle
NB.   iA   - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite, the inversion of A
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite, represented in the factored form by L
NB.          or U
NB.   n    ≥ 0, the size of A and iA
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dpotri=: 4 : 0
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  n=. # y
  cdrc=. dpotri_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)

zpotri=: 4 : 0
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  cdrc=. zpotri_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)
