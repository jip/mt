require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite tridiagonal matrix using
NB.   the factorization computed by xPTTRF
NB.
NB. Syntax:
NB.   X=.      dpttrs dD ; eT1 ; B
NB.   X=. uplo zpttrs dD ; eT1 ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.           specifies what is eT1:
NB.            'L' - eT1 is the subdiagonal of L1, the form
NB.                  is:
NB.                    L1 * D * L1^H = A
NB.            'U' - eT1 is the superdiagonal of U1, the form
NB.                  is:
NB.                    U1^H * D * U1 = A
NB.   dD   - n-vector, real, diagonal of D
NB.   eT1  - (n-1)-vector, subdiagonal of L1, or the
NB.          superdiagonal of U1
NB.   B    - n×nrhs-matrix, RHS
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite tridiagonal represented in factored
NB.          form by dD and eT1
NB.   D    - n×n-matrix, diagonal, with dD as diagonal
NB.          elements
NB.   L1   - n×n-matrix, the unit lower bidiagonal
NB.   U1   - n×n-matrix, the unit upper bidiagonal
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - no check for positive definiteness

dpttrs=: 3 : 0
  'd e B'=. y
  'n nrhs'=. $ B
  assert. (isvector_jlapack2_ , isreal_jlapack2_ ,     n  = #) d
  assert. (isvector_jlapack2_ , isreal_jlapack2_ , (<: n) = #) e
  assert. (ismatrix_jlapack2_ , isreal_jlapack2_             ) B
  cdrc=. dpttrs_jlapack2_ (, n) ; (, nrhs) ; d ; e ; (|: B) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 5 {:: cdrc
)

zpttrs=: 4 : 0
  'd e B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (isvector_jlapack2_ , isreal_jlapack2_ ,     n  = #) d
  assert. (isvector_jlapack2_ ,                    (<: n) = #) e
  assert.  ismatrix_jlapack2_                                  B
  cdrc=. zpttrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; d ; e ; (|: B) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  |: 6 {:: cdrc
)
