require 'math/lapack2'

NB. Description:
NB.   Compute the factorization:
NB.     L1 * D * L1^H = A
NB.   the same as
NB.     U1^H * D * U1 = A
NB.   of a Hermitian (symmetric) positive definite
NB.   tridiagonal matrix
NB.
NB. Syntax:
NB.   'dD eT1'=. xpttrf dA ; eA
NB. where
NB.   dA  - n-vector, real, diagonal of A
NB.   eA  - (n-1)-vector, subdiagonal of A
NB.   dD  - n-vector, real, diagonal of D
NB.   eT1 - (n-1)-vector, subdiagonal of L1, the same as
NB.         superdiagonal of U1
NB.   A   - n×n-matrix, the Hermitian (symmetric) positive
NB.         definite tridiagonal
NB.   D   - n×n-matrix, diagonal, with dD as diagonal
NB.         elements
NB.   L1  - n×n-matrix, the unit lower bidiagonal
NB.   U1  - n×n-matrix, the unit upper bidiagonal
NB.   n   ≥ 0, the size of A, D, L1 and U1
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - no check for positive definiteness

dpttrf=: 3 : 0
  'd e'=. y
  n=. # d
  assert.  isvector_jlapack2_               d
  assert. (isvector_jlapack2_ , (<: n) = #) e
  cdrc=. dpttrf_jlapack2_ (, n) ; d ; e ; , _1
  assert. 0 = _1 {:: cdrc
  2 3 { cdrc
)

zpttrf=: 3 : 0
  'd e'=. y
  n=. # d
  assert.  isvector_jlapack2_               d
  assert. (isvector_jlapack2_ , (<: n) = #) e
  cdrc=. zpttrf_jlapack2_ (, n) ; d ; e ; , _1
  assert. 0 = _1 {:: cdrc
  2 3 { cdrc
)
