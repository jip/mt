require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite tridiagonal matrix, and
NB.   factor the last one as:
NB.     L1 * D * L1^H = A
NB.
NB. Syntax:
NB.   'dD eT1 X'=. xptsv dA ; eA ; B
NB. where
NB.   dA  - n-vector, real, diagonal of A
NB.   eA  - (n-1)-vector, subdiagonal of A
NB.   B    - n×nrhs-matrix, RHS
NB.   dD   - n-vector, real, diagonal of D
NB.   eT1  - (n-1)-vector, subdiagonal of L1, or the
NB.          superdiagonal of U1
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite tridiagonal represented by dA and eA
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

dptsv=: 3 : 0
  'd e B'=. y
  'n nrhs'=. $ B
  assert. (isvector_jlapack2_ *. isreal_jlapack2_ *.     n  = #) d
  assert. (isvector_jlapack2_ *. isreal_jlapack2_ *. (<: n) = #) e
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_              ) B
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
  cdrc=. dptsv_jlapack2_ (, n) ; (, nrhs) ; d ; e ; (|: B) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 4 5 { cdrc  NB. (|:) doesn't affect to d and e
)

zptsv=: 3 : 0
  'd e B'=. y
  'n nrhs'=. $ B
  assert. (isvector_jlapack2_ *. isreal_jlapack2_ *.     n  = #) d
  assert. (isvector_jlapack2_ *.                     (<: n) = #) e
  assert.  ismatrix_jlapack2_                                 B
  select. 3!:0 d
    case. JCMPX do. d=. 9 o. d
    case. JFL   do.
    case.       do. d=. d + 0.0
  end.
  if. JCMPX ~: 3!:0 e do. e=. e + 0j0 end.
  cdrc=. zptsv_jlapack2_ (, n) ; (, nrhs) ; d ; e ; (|: B) ; (, 1 >. n) ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 4 5 { cdrc  NB. (|:) doesn't affect to d and e
)
