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
  assert. (isvector_jlapack2_ *. isreal_jlapack2_              ) d
  assert. (isvector_jlapack2_ *. isreal_jlapack2_ *. (<: n) = #) e
  select. 3!:0 d
    case. JCMPX do. d=. 9 o. d
    case. JFL   do.
    case.       do. d=. d + 0.0
  end.
  select. 3!:0 e
    case. JCMPX do. e=. 9 o. e
    case. JFL   do.
    case.       do. e=. e + 0.0
  end.
  cdrc=. dpttrf_jlapack2_ (, n) ; d ; e ; , _1
  assert. 0 = _1 {:: cdrc
  2 3 { cdrc
)

zpttrf=: 3 : 0
  'd e'=. y
  n=. # d
  assert. (isvector_jlapack2_ *. isreal_jlapack2_) d
  assert. (isvector_jlapack2_ *. (<: n) = #      ) e
  select. 3!:0 d
    case. JCMPX do. d=. 9 o. d
    case. JFL   do.
    case.       do. d=. d + 0.0
  end.
  if. JCMPX ~: 3!:0 e do. e=. e + 0j0 end.
  cdrc=. zpttrf_jlapack2_ (, n) ; d ; e ; , _1
  assert. 0 = _1 {:: cdrc
  2 3 { cdrc
)
